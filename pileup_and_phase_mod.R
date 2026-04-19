library(glue)
library(stringr)
library(data.table)
library(dplyr)
library(vcfR)
library(Matrix)
library(numbat)
library(peakRAM)
library(GenomicRanges)
library(R.utils)
library(optparse)

parser = OptionParser(description='Run input generation step using batching strategy')
parser = add_option(parser, '--sampleID', default = 'subject', type = "character", help = "Individual label. One per run.")
args <- parse_args(parser)

sampleID <- args$sampleID
ori_outdir <- paste0("/ix1/alee/LO_LAB/Personal/Rick/Haoyu/numbat/simulation_tm/preprocess/", sampleID, "/out")
gmap="/ix1/alee/LO_LAB/Personal/Rick/clone_review/Reference_genome/numbat_ref/genetic_map_hg38_withX.txt.gz"
genome = ifelse(str_detect(gmap, 'hg19'), 'hg19', 'hg38')
message(paste0('Using genome version: ', genome))
samples = str_split(sampleID, ',')[[1]]
# label here can only be one sample!! Not same as original Numbat (for testing only)
label = sampleID

# where to write the timing/memory file
out_file_dir <- paste0("/ix1/alee/LO_LAB/Personal/Rick/Haoyu/numbat/simulation_tm/preprocess_mod/", sampleID, "/")
dir.create(out_file_dir, showWarnings = FALSE, recursive = TRUE)
setwd(out_file_dir)
out_file <- glue("{out_file_dir}/timing_summary.txt")

#-------------------------------------------------------------------------------------
#' Preprocess allele data
#'
#' @param gt_col character Sample label
#' @param vcf_pu dataframe Pileup VCF from cell-snp-lite
#' @param vcf_phased dataframe Phased VCF from eagle2
#' @param AD dgTMatrix Alt allele depth matrix from pileup
#' @param DP dgTMatrix Total allele depth matrix from pileup
#' @param barcodes vector List of barcodes from pileup
#' @param gtf dataframe Transcript GTF
#' @param gmap dataframe Genetic map
#' @return dataframe Allele counts by cell
#' @keywords internal
preprocess_allele_chr = function(
    gt_col,
    vcf_pu,
    vcf_phased,
    AD,
    DP,
    barcodes,
    gtf,
    gmap
) {
  
  # pileup VCF
  vcf_pu = vcf_pu %>%
    mutate(INFO = str_remove_all(INFO, '[:alpha:]|=')) %>%
    tidyr::separate(col = 'INFO', into = c('AD', 'DP', 'OTH'), sep = ';') %>%
    mutate_at(c('AD', 'DP', 'OTH'), as.integer) %>%
    mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_'))
  
  # pileup count matrices
  DP = as.data.frame(Matrix::summary(DP)) %>%
    mutate(
      cell = barcodes[j],
      snp_id = vcf_pu$snp_id[i]
    ) %>%
    select(-i, -j) %>%
    dplyr::rename(DP = x) %>%
    select(cell, snp_id, DP)
  data.table::setDT(DP)
  
  
  AD = as.data.frame(Matrix::summary(AD)) %>%
    mutate(
      cell = barcodes[j],
      snp_id = vcf_pu$snp_id[i]
    ) %>%
    select(-i, -j) %>%
    dplyr::rename(AD = x) %>%
    select(cell, snp_id, AD)
  data.table::setDT(AD)
  
  # join AD+DP
  setkey(DP, cell, snp_id)
  setkey(AD, cell, snp_id)
  df = AD[DP]
  df[is.na(AD), AD := 0L]
  
  # join vcf_pu
  vcf_pu_dt = as.data.table(vcf_pu)
  vcf_pu_dt[, `:=`(AD_all = AD, DP_all = DP, OTH_all = OTH)]
  vcf_pu_dt = vcf_pu_dt[, .(snp_id, CHROM, POS, REF, ALT, AD_all, DP_all, OTH_all)]
  
  setkey(vcf_pu_dt, snp_id)
  setkey(df, snp_id)
  df = vcf_pu_dt[df]
  
  df[, `:=`(
    AR = AD / DP,
    AR_all = AD_all / DP_all
  )]
  df = df[DP_all > 1 & OTH_all == 0]
  
  # vcf has duplicated records sometimes
  df = unique(df, by = c("cell", "snp_id"))
  
  # df = DP %>% left_join(AD, by = c("cell", "snp_id")) %>%
  #   mutate(AD = ifelse(is.na(AD), 0, AD))
  #
  # df = df %>% left_join(
  #   vcf_pu %>% dplyr::rename(AD_all = AD, DP_all = DP, OTH_all = OTH),
  #   by = 'snp_id')
  #
  # df = df %>% mutate(
  #   AR = AD/DP,
  #   AR_all = AD_all/DP_all
  # )
  #
  # df = df %>% dplyr::filter(DP_all > 1 & OTH_all == 0)
  #
  # # vcf has duplicated records sometimes
  # df = df %>% distinct()
  #
  # df = df %>% mutate(
  #   snp_index = as.integer(factor(snp_id, unique(snp_id))),
  #   cell_index = as.integer(factor(cell, base::sample(unique(cell))))
  # )
  
  # phased VCF
  vcf_phased = as.data.table(vcf_phased)
  vcf_phased[, snp_id := paste(CHROM, POS, REF, ALT, sep = '_')]
  vcf_phased[, GT := .SD[[gt_col]], .SDcols = gt_col]
  #
  # vcf_phased = vcf_phased %>%
  #   mutate(snp_id = paste(CHROM, POS, REF, ALT, sep = '_')) %>%
  #   mutate(GT = get(sample))
  
  # annotate SNP by gene
  hits = GenomicRanges::findOverlaps(
    GenomicRanges::GRanges(seqnames = vcf_phased$CHROM,
                           ranges = IRanges::IRanges(start = vcf_phased$POS, end = vcf_phased$POS)),
    GenomicRanges::GRanges(seqnames = gtf$CHROM,
                           ranges = IRanges::IRanges(start = gtf$gene_start, end = gtf$gene_end))
  )
  
  if (length(hits) > 0) {
    olap = as.data.table(as.data.frame(hits))
    setnames(olap, c("snp_index", "gene_index"))
    vcf_phased[, snp_index := .I]
    gtf = as.data.table(gtf)
    gtf[, gene_index := .I]
    
    olap = merge(olap, vcf_phased[, .(snp_index, snp_id)], by = "snp_index")
    olap = merge(olap, gtf[, .(gene_index, gene, gene_start, gene_end)], by = "gene_index")
    setorder(olap, snp_index, gene)
    olap = unique(olap, by = "snp_id")
    olap = olap[, .(snp_id, gene, gene_start, gene_end)]
  } else {
    olap = data.table(snp_id = character(), gene = character(),
                      gene_start = integer(), gene_end = integer())
  }
  
  setkey(vcf_phased, snp_id)
  setkey(olap, snp_id)
  vcf_phased = olap[vcf_phased]
  
  # overlap_transcript = GenomicRanges::findOverlaps(
  #   vcf_phased %>% {GenomicRanges::GRanges(
  #     seqnames = .$CHROM,
  #     IRanges::IRanges(start = .$POS,
  #                      end = .$POS)
  #   )},
  #   gtf %>% {GenomicRanges::GRanges(
  #     seqnames = .$CHROM,
  #     IRanges::IRanges(start = .$gene_start,
  #                      end = .$gene_end)
  #   )}
  # ) %>%
  #   as.data.frame() %>%
  #   setNames(c('snp_index', 'gene_index')) %>%
  #   left_join(
  #     vcf_phased %>% mutate(snp_index = 1:n()) %>%
  #       select(snp_index, snp_id),
  #     by = c('snp_index')
  #   ) %>%
  #   left_join(
  #     gtf %>% mutate(gene_index = 1:n()),
  #     by = c('gene_index')
  #   ) %>%
  #   arrange(snp_index, gene) %>%
  #   distinct(snp_index, `.keep_all` = TRUE)
  #
  # vcf_phased = vcf_phased %>%
  #   left_join(
  #     overlap_transcript %>% select(snp_id, gene, gene_start, gene_end),
  #     by = c('snp_id')
  #   )
  
  # annotate SNPs by genetic map
  map_hits = GenomicRanges::findOverlaps(
    GenomicRanges::GRanges(seqnames = vcf_phased$CHROM,
                           ranges = IRanges::IRanges(start = vcf_phased$POS, end = vcf_phased$POS)),
    GenomicRanges::GRanges(seqnames = gmap$CHROM,
                           ranges = IRanges::IRanges(start = gmap$start, end = gmap$end))
  )
  
  if (length(map_hits) > 0) {
    map_tbl = as.data.table(as.data.frame(map_hits))
    setnames(map_tbl, c("marker_index", "map_index"))
    
    vcf_phased[, marker_index := .I]
    gmap = as.data.table(gmap)
    gmap[, map_index := .I]
    
    marker_map = merge(map_tbl, vcf_phased[, .(marker_index, snp_id)], by = "marker_index")
    marker_map = merge(marker_map, gmap[, .(map_index, start, cM)], by = "map_index", all.x = TRUE)
    
    setorder(marker_map, marker_index, -start)
    marker_map = unique(marker_map, by = "marker_index")
    marker_map = marker_map[, .(snp_id, cM)]
  } else {
    marker_map = data.table(snp_id = character(), cM = numeric())
  }
  
  setkey(vcf_phased, snp_id)
  setkey(marker_map, snp_id)
  vcf_phased = marker_map[vcf_phased]
  
  # marker_map = GenomicRanges::findOverlaps(
  #   vcf_phased %>% {GenomicRanges::GRanges(
  #     seqnames = .$CHROM,
  #     IRanges::IRanges(start = .$POS,
  #                      end = .$POS)
  #   )},
  #   gmap %>% {GenomicRanges::GRanges(
  #     seqnames = .$CHROM,
  #     IRanges::IRanges(start = .$start,
  #                      end = .$end)
  #   )}
  # ) %>%
  #   as.data.frame() %>%
  #   setNames(c('marker_index', 'map_index')) %>%
  #   left_join(
  #     vcf_phased %>% mutate(marker_index = 1:n()) %>%
  #       select(marker_index, snp_id),
  #     by = c('marker_index')
  #   ) %>%
  #   left_join(
  #     gmap %>% mutate(map_index = 1:n()),
  #     by = c('map_index')
  #   ) %>%
  #   arrange(marker_index, -start) %>%
  #   distinct(marker_index, `.keep_all` = TRUE) %>%
  #   select(snp_id, cM)
  #
  # vcf_phased = vcf_phased %>%
  #   left_join(marker_map, by = 'snp_id')
  
  # add annotation to cell counts
  setkey(df, snp_id)
  df = vcf_phased[, .(snp_id, gene, GT, cM)][df]
  
  df = df[GT %in% c('1|0', '0|1'),
          .(cell, snp_id, CHROM, POS, cM, REF, ALT, AD, DP, GT, gene)]
  
  df[, CHROM := factor(CHROM, levels = unique(vcf_pu$CHROM))]
  
  # df = df %>% left_join(vcf_phased %>% select(snp_id, gene, GT, cM), by = 'snp_id')
  # df = df %>% mutate(CHROM = factor(CHROM, unique(CHROM)))
  # df = df %>% select(cell, snp_id, CHROM, POS, cM, REF, ALT, AD, DP, GT, gene)
  #
  # # only keep hets
  # df = df %>% filter(GT %in% c('1|0', '0|1'))
  
  return(df)
}
#-------------------------------------------------------------------------------------

# run function
res <- peakRAM({
  
  if (genome == 'hg19') {
    gtf = gtf_hg19
  } else {
    gtf = gtf_hg38
  }
  
  genetic_map = fread(gmap) %>%
    setNames(c('CHROM', 'POS', 'rate', 'cM')) %>%
    group_by(CHROM) %>%
    mutate(
      start = POS,
      end = c(POS[2:length(POS)], POS[length(POS)])
    ) %>%
    ungroup()
  
  for (sample in samples) {
    
    outfile <- glue("{out_file_dir}/{sample}_allele_counts.tsv")
    fwrite(
      data.table(cell=character(), snp_id=character(), CHROM=character(),
                 POS=integer(), cM=numeric(), REF=character(), ALT=character(),
                 AD=integer(), DP=integer(), GT=character(), gene=character()),
      outfile, sep = '\t', append = FALSE
    )
    
    pu_dir = glue('{ori_outdir}/pileup/{sample}')
    
    # pileup VCF
    vcf_pu <- fread(glue('{pu_dir}/cellSNP.base.vcf'), skip = '#CHROM')
    setnames(vcf_pu, old = "#CHROM", new = "CHROM")
    vcf_pu[, CHROM := sub("^chr", "", CHROM)]
    
    
    # count matrices (SNP x cell)
    # each nonzero entry is represented explicitly by a triplet (i, j, x)
    # i-th SNP in j-th cell has x count.
    AD = readMM(glue('{pu_dir}/cellSNP.tag.AD.mtx'))
    DP = readMM(glue('{pu_dir}/cellSNP.tag.DP.mtx'))
    
    
    cell_barcodes = fread(glue('{pu_dir}/cellSNP.samples.tsv'), header = F) %>% pull(V1)
    
    # chunk by chromosome
    chroms <- intersect(as.character(1:22), unique(vcf_pu$CHROM))
    
    for (chr in chroms) {
      message(glue('Processing chr {chr} ...'))
      
      # chunk pileup VCF by chr
      vcf_pu_chr <-  vcf_pu %>% filter(CHROM == chr)
      idx_chr <- which(vcf_pu$CHROM == chr)
      
      # chunk AD, DP by chr
      AD_chr <- AD[idx_chr, , drop = FALSE]
      DP_chr <- DP[idx_chr, , drop = FALSE]
      # AD_chr <- as(AD_chr, "dgTMatrix")
      # DP_chr <- as(DP_chr, "dgTMatrix")
      
      # chunk gtf, genetic map by chr
      gtf_chr <- gtf %>% filter(CHROM == chr)
      gmap_chr <- genetic_map %>% filter(CHROM == chr)
      
      # read in phased VCF per chr
      vcf_file = glue('{ori_outdir}/phasing/{label}_chr{chr}.phased.vcf.gz')
      if (file.exists(vcf_file)) {
        vcf_phased_chr <- fread(vcf_file, skip = '#CHROM')
        setnames(vcf_phased_chr, old = "#CHROM", new = "CHROM")
        vcf_phased_chr[, CHROM := sub("^chr", "", CHROM)]
        
        # vcf_phased_chr <- fread(vcf_file, skip = '#CHROM') %>%
        #   rename(CHROM = `#CHROM`) %>%
        #   mutate(CHROM = str_remove(CHROM, 'chr'))
      } else {
        stop('Phased VCF not found')
      }
      
      
      df_chr = preprocess_allele_chr(
        gt_col = label,
        vcf_pu = vcf_pu_chr,
        vcf_phased = vcf_phased_chr,
        AD = AD_chr,
        DP = DP_chr,
        barcodes = cell_barcodes,
        gtf = gtf_chr,
        gmap = gmap_chr
      ) %>%
        filter(GT %in% c('1|0', '0|1'))
      
      # append results
      fwrite(df_chr, outfile, sep = '\t', append = TRUE)
      rm(vcf_pu_chr, vcf_phased_chr, AD_chr, DP_chr, gtf_chr, gmap_chr, df_chr)
      gc()
    }
    
    R.utils::gzip(outfile, overwrite = TRUE)
    file.remove(outfile)
    
    
    # vcf_phased = lapply(1:22, function(chr) {
    #   vcf_file = glue('{outdir}/phasing/{label}_chr{chr}.phased.vcf.gz')
    #   if (file.exists(vcf_file)) {
    #     fread(vcf_file, skip = '#CHROM') %>%
    #       rename(CHROM = `#CHROM`) %>%
    #       mutate(CHROM = str_remove(CHROM, 'chr'))
    #   } else {
    #     stop('Phased VCF not found')
    #   }
    # }) %>%
    #   Reduce(rbind, .) %>%
    #   mutate(CHROM = factor(CHROM, unique(CHROM)))
    #
    #     df = preprocess_allele_mod(
    #       sample = label,
    #       vcf_pu = vcf_pu,
    #       vcf_phased = vcf_phased,
    #       AD = AD,
    #       DP = DP,
    #       barcodes = cell_barcodes,
    #       gtf = gtf,
    #       gmap = genetic_map
    #     ) %>%
    #       filter(GT %in% c('1|0', '0|1'))
    #
    #     fwrite(df, glue('{out_file_dir}/{sample}_allele_counts.tsv.gz'), sep = '\t')
    
  }
  
})


# write out timing and memory result
write.table(
  res[, c("Elapsed_Time_sec", "Peak_RAM_Used_MiB")],
  file = out_file,
  sep = "\t", quote = FALSE, row.names = FALSE,
  col.names = c("elapsed_sec", "peak_ram_mib")
)
