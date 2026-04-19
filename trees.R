# from phangorn
#' UPGMA and WPGMA clustering
#'
#' @param D A distance matrix.
#' @param method The agglomeration method to be used. This should be (an
#' unambiguous abbreviation of) one of "ward", "single", "complete", "average",
#' "mcquitty", "median" or "centroid". The default is "average".
#' @param \dots Further arguments passed to or from other methods.
upgma <- function(D, method = "average", ...) {
  DD <- as.dist(D)
  hc <- hclust(DD, method = method, ...)
  result <- ape::as.phylo(hc)
  result <- reorder(result, "postorder")
  result
}

#' Mark the tumor lineage of a phylogeny
#' @param gtree tbl_graph Single-cell phylogeny
#' @return tbl_graph Phylogeny annotated with tumor versus normal compartment
#' @keywords internal
mark_tumor_lineage = function(gtree) {

  mut_nodes = gtree %>%
    activate(nodes) %>%
    filter(!is.na(site)) %>%
    as.data.frame() %>%
    pull(id)

  mut_burdens = lapply(
    mut_nodes,
    function(node) {
      gtree %>%
        activate(nodes) %>%
        mutate(
          mut_burden = ifelse(GT == '', 0, str_count(GT, ',') + 1)
        ) %>%
        ungroup() %>%
        mutate(seq = bfs_rank(root = node)) %>%
        data.frame %>%
        filter(leaf & seq > 0) %>%
        pull(mut_burden) %>%
        sum
    }
  )

  tumor_root = mut_nodes[which.max(mut_burdens)]

  gtree = gtree %>%
    activate(nodes) %>%
    mutate(
      seq = bfs_rank(root = tumor_root),
      compartment = ifelse(seq > 0, 'tumor', 'normal'),
      is_tumor_root = tumor_root == id
    )

  compartment_dict = gtree %>%
    activate(nodes) %>%
    as.data.frame() %>%
    {setNames(.$compartment, .$id)}

  gtree = gtree %>%
    activate(edges) %>%
    mutate(compartment = compartment_dict[to])

  return(gtree)

}


#' Find maximum lilkelihood assignment of mutations on a tree
#' @param tree phylo Single-cell phylogenetic tree
#' @param P matrix Genotype probability matrix
#' @return list Mutation
#' @keywords internal
get_tree_post = function(tree, P) {

  sites = colnames(P)
  n = nrow(P)
  tree_stats = score_tree(tree, P, get_l_matrix = TRUE)

  # l_matrix to find the place best put the segments.

  l_matrix = as.data.frame(tree_stats$l_matrix)

  colnames(l_matrix) = sites
  rownames(l_matrix) = c(tree$tip.label, paste0('Node', 1:tree$Nnode))

  # place the mutations to the its best place on phylogeny tree.
  gtree = annotate_tree(tree, P)

  return(list('gtree' = gtree, 'l_matrix' = l_matrix))
}

#' Get a tidygraph tree with simplified mutational history.
#' @description Specify either max_cost or n_cut.
#' max_cost works similarly as h and n_cut works similarly as k in stats::cutree.
#' The top-level normal diploid clone is always included.
#' @param tree phylo Single-cell phylogenetic tree
#' @param P matrix Genotype probability matrix
#' @param max_cost numeric Likelihood threshold to collapse internal branches
#' @param n_cut integer Number of cuts on the phylogeny to define subclones
#' @return tbl_graph Phylogeny annotated with branch lengths and mutation events
#' @export
get_gtree = function(tree, P, n_cut = 0, max_cost = 0) {

  # calculate L matrix and create tidygraph object
  # tree_post is still a phylogeny tree but with mutation annotated (intermediate containter)

  # ---------------------------------------
  # here, the tree should be built on centroids
  # turn P (single cell x site matrix) to P2 (centroids x sites matrix)
  # ---------------------------------------

  tree_post = get_tree_post(tree, P)


  # simplify mutational history
  # get_mut_graph: Convert a single-cell phylogeny with mutation placements into a mutation graph

  # G_m is not a cell‐tree, but rather a graph of genotype states (clones), where:
  # Vertices = distinct genotype nodes (each carrying one or more mutation “labels” and a cumulative GT string),
  # Edges = “mutation happened here → next mutation happened there” relationships.
  #
  # (root)
  #   |
  #   [M1]
  #   |
  #   (n2)
  # /   \
  # [M2]  [M3]
  # /       \
  # (n3)     (n4)
  #
  # Vertices:
  #   | id | label  | GT        | clone |
  #   | -- | ------ | --------- | ----- |
  #   | 1  | `""`   | `""`      | 1     |
  #   | 2  | `"M1"` | `"M1"`    | 2     |
  #   | 3  | `"M2"` | `"M1,M2"` | 3     |
  #   | 4  | `"M3"` | `"M1,M3"` | 4     |
  #
  #
  # Edges:
  #   | from | to |
  #   | ---- | -- |
  #   | 1    | 2  |
  #   | 2    | 3  |
  #   | 2    | 4  |
  #
  # simplify_history: prune the tree
  # Weakly supported mutation nodes (low penalty to reassign) get merged with their better‐supported neighbor.
  # End up with a smaller graph of surviving mutation nodes, each with a label that may now represent multiple originally distinct events.
  # That simplified graph is then fed to label_genotype(), which recomputes for each merged vertex the cumulative genotype string (GT) and a new clone ID.
  #

  gtree <- tree_post$gtree

  G_m = get_mut_graph(gtree)
  G_m = simplify_history(G_m, tree_post$l_matrix, max_cost = max_cost, n_cut = n_cut)
  G_m = label_genotype(G_m)

  mut_nodes = G_m %>% igraph::as_data_frame('vertices') %>%
    select(name = node, site = label, clone = clone, GT = GT)

  # update tree
  # mut_to_tree: prune the phylogeny tree (gtree) according to the pruned mutation graph.
  gtree = mut_to_tree(gtree, mut_nodes)

  # gtree <- gtree %>%
  #   activate(edges) %>%
  #   mutate(.tidygraph_edge_index = row_number()) %>%
  #   activate(nodes)

  # mark_tumor_lineage:
  # dentify all mutation‐event nodes (!is.na(site)).
  # For each one, sum the total number of mutation calls downstream of it (i.e. on all descendant leaves), giving each node a “mutation burden.”
  # Choose the node with the largest burden → call this the tumor root.
  # Run a breadth‐first search from that tumor rooT.

  gtree = mark_tumor_lineage(gtree)

  # gtree <- gtree %>%
  #   activate(edges) %>%
  #   mutate(.tidygraph_edge_index = row_number()) %>%
  #   activate(nodes)


  # NOTICE: gtree is a clone tree
  # (a phylogenetic tree whose vertices represent the final set of clones (genotypes), not each individual cell):
  #   |  id |  site  | clone |     GT    |  leaf | compartment | is\_tumor\_root | seq |
  #   | :-: | :----: | :---: | :-------: | :---: | :---------- | :-------------- | :-: |
  #   |  1  |   NA   |   1   |    `""`   | FALSE | normal      | FALSE           |  0  |
  #   |  3  | `"M1"` |   2   |   `"M1"`  | FALSE | tumor       | TRUE            |  0  |
  #   |  5  | `"M2"` |   3   | `"M1,M2"` |  TRUE | tumor       | FALSE           |  1  |
  # seq is the BFS rank


  return(gtree)

}

#' Annotate the direct upstream or downstream node on the edges
#'
#' @param G igraph Mutation graph
#' @return igraph Mutation graph
#' @keywords internal
transfer_links = function(G) {

  edge_df = G %>% igraph::as_data_frame('edges') %>%
    left_join(
      G %>% igraph::as_data_frame('vertices') %>% select(from_node = node, id),
      by = c('from' = 'id')
    ) %>%
    left_join(
      G %>% igraph::as_data_frame('vertices') %>% select(to_node = node, id),
      by = c('to' = 'id')
    )

  E(G)$from_node = edge_df$from_node
  E(G)$to_node = edge_df$to_node

  return(G)
}

#' Label the genotypes on a mutation graph
#'
#' @param G igraph Mutation graph
#' @return igraph Mutation graph
#' @keywords internal
label_genotype = function(G) {

  id_to_label = igraph::as_data_frame(G, 'vertices') %>% {setNames(.$label, .$id)}

  # for some reason, the output from all_simple_path is out of order if supplied directly
  # V(G)$GT = igraph::all_simple_paths(G, from = 1) %>%
  V(G)$GT = lapply(
    2:length(V(G)),
    function(v) {dplyr::first(igraph::all_simple_paths(G, from = 1, to = v), default = NULL)}
  ) %>%
    purrr::map(as.character) %>%
    purrr::map(function(x) {
      muts = id_to_label[x]
      muts = muts[muts != '']
      paste0(muts, collapse = ',')
    }) %>%
    c(id_to_label[[1]],.) %>%
    as.character

  visit_order = setNames(1:length(V(G)), as.numeric(igraph::dfs(G, root = 1)$order))
  V(G)$clone = visit_order[as.character(as.numeric(V(G)))]

  return(G)
}

#' Annotate the direct upstream or downstream mutations on the edges
#'
#' @param G igraph Mutation graph
#' @return igraph Mutation graph
#' @keywords internal
label_edges = function(G) {

  edge_df = G %>% igraph::as_data_frame('edges') %>%
    left_join(
      G %>% igraph::as_data_frame('vertices') %>% select(from_label = label, id),
      by = c('from' = 'id')
    ) %>%
    left_join(
      G %>% igraph::as_data_frame('vertices') %>% select(to_label = label, id),
      by = c('to' = 'id')
    ) %>%
    mutate(label = paste0(from_label, '->', to_label))

  E(G)$label = edge_df$label
  E(G)$from_label = edge_df$from_label
  E(G)$to_label = edge_df$to_label

  return(G)
}

#' Merge adjacent set of nodes
#'
#' @param G igraph Mutation graph
#' @param vset vector Set of adjacent vertices to merge
#' @return igraph Mutation graph
#' @keywords internal
contract_nodes = function(G, vset, node_tar = NULL, debug = FALSE) {

  vset = unlist(vset)

  if (length(vset) == 1) {
    return(G)
  }

  # reorder the nodes according to graph
  vorder = V(G)$label[igraph::dfs(G, root = 1)$order]
  vset = vorder[vorder %in% vset]

  vset_ids = V(G)[label %in% vset]

  ids_new = 1:vcount(G)

  # the indices before do not change
  ids_new[vset_ids] = min(vset_ids)
  # indices after might need to be reset
  if (max(vset_ids) != vcount(G)) {
    ids_new[(max(vset_ids)+1):length(ids_new)] = ids_new[(max(vset_ids)+1):length(ids_new)] - length(vset_ids) + 1
  }

  G = G %>% igraph::contract(
    ids_new,
    vertex.attr.comb = list(label = function(x){paste0(sort(x), collapse = ',')}, node = "first", "ignore")
  )

  if (!is.null(node_tar)) {
    V(G)[min(vset_ids)]$node = node_tar
  }

  V(G)$id = 1:vcount(G)

  G = igraph::simplify(G)

  if (debug) {
    return(G)
  }

  G = label_edges(G)

  return(G)

}

#' Simplify the mutational history based on likelihood evidence
#'
#' @param G igraph Mutation graph
#' @param l_matrix matrix Mutation placement likelihood matrix (node by mutation)
#' @return igraph Mutation graph
#' @keywords internal
simplify_history = function(G, l_matrix, max_cost = 150, n_cut = 0, verbose = TRUE) {

  if (n_cut > 0) {
    max_cost = Inf
  }

  for (i in 1:ecount(G)) {

    move_opt = get_move_opt(G, l_matrix)

    if (move_opt$cost < max_cost & ecount(G) > n_cut-1) {
      log_message("Start collapsing edges...")

      if (move_opt$direction == 'up') {
        G = G %>% contract_nodes(c(move_opt$from_label, move_opt$to_label), move_opt$from_node) %>% transfer_links()
        msg = glue('opt_move:{move_opt$to_label}->{move_opt$from_label}, cost={signif(move_opt$cost,3)}')
      } else {
        G = G %>% contract_nodes(c(move_opt$from_label, move_opt$to_label), move_opt$to_node) %>% transfer_links()
        msg = glue('opt_move:{move_opt$from_label}->{move_opt$to_label}, cost={signif(move_opt$cost,3)}')
      }

      log_info(msg)

    } else {
      break()
    }
  }

  return(G)
}

#' Get the cost of a mutation reassignment
#'
#' @param muts character Mutations dlimited by comma
#' @param node_ori character Name of the "from" node
#' @param node_tar character Name of the "to" node
#' @return numeric Likelihood cost of the mutation reassignment
#' @keywords internal
get_move_cost = function(muts, node_ori, node_tar, l_matrix) {

  if (muts == '') {
    return(Inf)
  }

  if (str_detect(muts, ',')) {
    muts = unlist(str_split(muts, ','))
  }

  sum(l_matrix[node_ori, muts] - l_matrix[node_tar, muts])
}

#' Get the least costly mutation reassignment
#'
#' @param G igraph Mutation graph
#' @param l_matrix matrix Likelihood matrix of mutation placements
#' @return numeric Lieklihood cost of performing the mutation move
#' @keywords internal
get_move_opt = function(G, l_matrix) {

  move_opt = G %>% igraph::as_data_frame('edges') %>%
    group_by(from) %>%
    mutate(n_sibling = n()) %>%
    ungroup() %>%
    rowwise() %>%
    mutate(
      up = get_move_cost(to_label, to_node, from_node, l_matrix),
      down = get_move_cost(from_label, from_node, to_node, l_matrix)
    ) %>%
    ungroup() %>%
    # prevent a down move if branching. Technically it's fine but graph has to be modified correctly
    mutate(down = ifelse(n_sibling > 1, Inf, down)) %>%
    as.data.table %>%
    data.table::melt(measure.vars = c('up', 'down'), variable.name = 'direction', value.name = 'cost') %>%
    arrange(cost) %>%
    head(1)

  return(move_opt)
}

#' Get ordered tips from a tree
#' @keywords internal
get_ordered_tips = function(tree) {
  is_tip <- tree$edge[,2] <= length(tree$tip.label)
  ordered_tips <- tree$edge[is_tip, 2]
  tree$tip.label[ordered_tips]
}



# # Unique root in 'phylo' numbering (tips = 1..n, internals = n+1..n+Nnode)
# phylo_root_id <- function(tree) {
#   parents  <- tree$edge[,1]
#   children <- tree$edge[,2]
#   r <- setdiff(parents, children)
#   if (length(r) != 1) stop("Cannot determine a unique root node.")
#   r[[1]]
# }
#
# # Strict postorder of internal nodes (> n) for upward DP
# postorder_internal <- function(tree) {
#   n <- length(tree$tip.label)
#   parents  <- tree$edge[,1]
#   children <- tree$edge[,2]
#   chlist <- split(children, parents)
#
#   out <- integer(0)
#   rec <- function(u) {
#     kids <- chlist[[as.character(u)]]
#     if (length(kids)) for (v in kids) rec(v)
#     out <<- c(out, u)
#   }
#
#   rec(phylo_root_id(tree))
#   out[out > n]
# }
#
#
# # score_tree for multifurcation tree
# score_tree_multifurcat <- function(tree, P, get_l_matrix = FALSE) {
#
#   tree <- reorder(tree, order = 'postorder')
#
#   n <- nrow(P)
#   m <- ncol(P)
#
#   logQ   <- matrix(nrow = tree$Nnode + n, ncol = m)
#   logP_0 <- log(1 - P)
#   logP_1 <- log(P)
#
#   logQ[1:n,] <- logP_1 - logP_0                # tips
#
#
#   children_dict <- allChildrenCPP(tree$edge)
#
#   # TRUE postorder of internal nodes (> n)
#   # may need to use faster way..
#   node_order <- postorder_internal(tree)
#
#   logQ <- CgetQ(logQ, children_dict, node_order)
#
#   if (get_l_matrix) {
#     l_matrix <- sweep(logQ, 2, colSums(logP_0), FUN = '+')
#     l_tree   <- sum(apply(l_matrix, 2, max))
#   } else {
#     l_matrix <- NULL
#     l_tree   <- sum(apply(logQ, 2, max)) + sum(logP_0)
#   }
#
#   list('l_tree' = l_tree, 'logQ' = logQ, 'l_matrix' = l_matrix)
# }
#
#
# # annotate_tree for multifurcation tree
# annotate_tree_multifurcat <- function(tree, P) {
#   sites <- colnames(P)
#   n     <- nrow(P)
#
#   tree_stats <- score_tree_multifurcat(tree, P, get_l_matrix = TRUE)
#   l_matrix   <- as.data.frame(tree_stats$l_matrix)
#   colnames(l_matrix) <- sites
#   rownames(l_matrix) <- c(tree$tip.label, paste0('Node', 1:tree$Nnode))
#
#   # mutation assignment on nodes
#   mut_nodes <- data.frame(
#     site       = sites,
#     node_phylo = apply(l_matrix, 2, which.max),
#     l          = apply(l_matrix, 2, max)
#   ) %>%
#     mutate(name = ifelse(node_phylo <= n, tree$tip.label[node_phylo], paste0('Node', node_phylo - n))) %>%
#     group_by(name) %>%
#     summarise(
#       site = paste0(sort(site), collapse = ','),
#       n_mut = n(),
#       l = sum(l),
#       .groups = 'drop'
#     )
#
#   #use phylo root id as character name for tidygraph bfs_*
#   root_phy <- phylo_root_id(tree)
#
#   gtree <- tree %>%
#     as_tbl_graph() %>%
#     mutate(
#       leaf  = node_is_leaf(),
#       # depth = bfs_dist(root = as.character(root_phy)),
#       # id    = as.integer(name)
#       root = node_is_root(),
#       depth = bfs_dist(root = 1),
#       id = row_number()
#     )
#
#   # leaf annotation for edges
#   gtree <- gtree %>%
#     activate(edges) %>%
#     select(-any_of(c('leaf'))) %>%
#     left_join(
#       gtree %>%
#         activate(nodes) %>%
#         data.frame() %>%
#         select(id, leaf),
#       by = c('to' = 'id')
#     )
#
#   # annotate with mutation placements
#   gtree <- mut_to_tree_multifurcat(gtree, mut_nodes)
#   log_message("safe here 5!")
#
#   gtree
# }
#
# # get_tree_post for multifurcation tree
# get_tree_post_multifurcat <- function(tree, P) {
#   sites <- colnames(P)
#   n     <- nrow(P)
#
#   tree_stats <- score_tree_multifurcat(tree, P, get_l_matrix = TRUE)
#   l_matrix   <- as.data.frame(tree_stats$l_matrix)
#
#   colnames(l_matrix) <- sites
#   rownames(l_matrix) <- c(tree$tip.label, paste0('Node', 1:tree$Nnode))
#   log_message("Safe here 2!")
#   gtree <- annotate_tree_multifurcat(tree, P)
#   return(list('gtree' = gtree, 'l_matrix' = l_matrix))
# }
#
#
# mark_tumor_lineage_multifurcat <- function(gtree) {
#
#   mut_nodes <- gtree %>%
#     activate(nodes) %>%
#     filter(!is.na(site)) %>%
#     as.data.frame() %>%
#     pull(id)
#
#   # Precompute per-leaf mutation burden once per candidate
#   mut_burdens <- lapply(
#     mut_nodes,
#     function(node_phy) {
#       gtree %>%
#         activate(nodes) %>%
#         mutate(mut_burden = ifelse(GT == '', 0, stringr::str_count(GT, ',') + 1)) %>%
#         ungroup() %>%
#         mutate(seq = bfs_rank(root = as.character(node_phy))) %>%
#         data.frame() %>%
#         dplyr::filter(leaf & seq > 0) %>%
#         dplyr::pull(mut_burden) %>%
#         sum()
#     }
#   )
#
#   tumor_root = mut_nodes[which.max(mut_burdens)]
#
#   gtree <- gtree %>%
#     activate(nodes) %>%
#     mutate(
#       seq         = bfs_rank(root = as.character(tumor_root)),
#       compartment = ifelse(seq > 0, 'tumor', 'normal'),
#       is_tumor_root = (id == tumor_root)
#     )
#
#   comp_dict <- gtree %>%
#     activate(nodes) %>%
#     as.data.frame() %>% { setNames(.$compartment, .$id) }
#
#   gtree %>%
#     activate(edges) %>%
#     mutate(compartment = comp_dict[to]) %>%
#     activate(nodes)
# }
#
#
# label_genotype_multifurcat <- function(G) {
#   vcnt <- igraph::vcount(G)
#   labels <- igraph::V(G)$label
#
#   # true root by in-degree 0
#   root_id <- which(igraph::degree(G, mode = "in") == 0)
#   if (length(root_id) != 1) stop("Mutation graph must have exactly one root.")
#
#   GT <- character(vcnt)
#   GT[root_id] <- labels[root_id]
#
#   for (v in setdiff(seq_len(vcnt), root_id)) {
#     p <- igraph::all_simple_paths(G, from = root_id, to = v, mode = "out")
#     if (length(p) == 0) next
#     path <- as.integer(p[[1]])
#     muts <- labels[path]
#     muts <- muts[!is.na(muts) & muts != ""]
#     GT[v] <- paste0(muts, collapse = ",")
#   }
#
#   igraph::V(G)$GT <- GT
#
#   # clone ids: DFS order from true root
#   dfs_order  <- igraph::dfs(G, root = root_id)$order
#   igraph::V(G)$clone <- match(seq_len(vcnt), dfs_order)
#
#   G
# }
#
#
#
# get_gtree_multifurcat <- function(tree, P, n_cut = 0, max_cost = 0) {
#
#   # calculate L matrix and create tidygraph object
#   # tree_post is still a phylogeny tree but with mutation annotated (intermediate containter)
#
#   log_message("Safe here!")
#   tree_post = get_tree_post_multifurcat(tree, P)
#   gtree <- tree_post$gtree
#
#   G_m = get_mut_graph(gtree)
#   G_m = simplify_history(G_m, tree_post$l_matrix, max_cost = max_cost, n_cut = n_cut)
#   G_m = label_genotype_multifurcat(G_m)
#
#   mut_nodes = G_m %>% igraph::as_data_frame('vertices') %>%
#     select(name = node, site = label, clone = clone, GT = GT)
#
#   # update tree
#   # mut_to_tree: prune the phylogeny tree (gtree) according to the pruned mutation graph.
#   gtree = mut_to_tree_multifurcat(gtree, mut_nodes)
#   gtree = mark_tumor_lineage_multifurcat(gtree)
#
#   return(gtree)
#
# }
#
#
# mut_to_tree_multifurcat <- function(gtree, mut_nodes) {
#
#   # transfer mutation to tree
#   gtree = gtree %>%
#     activate(nodes) %>%
#     select(-any_of(c('n_mut', 'l', 'site', 'clone'))) %>%
#     left_join(
#       mut_nodes %>%
#         mutate(n_mut = unlist(lapply(str_split(site, ','), length))) %>%
#         select(name, n_mut, site),
#       by = 'name'
#     ) %>%
#     mutate(n_mut = ifelse(is.na(n_mut), 0, n_mut))
#
#   # get branch length
#   gtree = gtree %>%
#     activate(edges) %>%
#     select(-any_of(c('length'))) %>%
#     left_join(
#       gtree %>%
#         activate(nodes) %>%
#         data.frame() %>%
#         select(id, length = n_mut),
#       by = c('to' = 'id')
#     ) %>%
#     mutate(length = ifelse(leaf, pmax(length, 0.2), length))
#
#   # label genotype on nodes
#   node_to_mut = gtree %>% activate(nodes) %>% data.frame() %>% {setNames(.$site, .$id)}
#
#   gtree = gtree %>%
#     activate(nodes) %>%
#     mutate(
#       .root = which(node_is_root())[1],
#       GT = vapply(
#         map_bfs(.root,
#                 .f = function(path, ...) {
#                   paste0(na.omit(node_to_mut[path$node]), collapse = ',')
#                 },
#                 .mode = "out"
#         ),
#         FUN = identity,
#         FUN.VALUE = character(1)
#       ),
#       last_mut = vapply(
#         map_bfs(.root,
#                 .f = function(path, ...) {
#                   past_muts = na.omit(node_to_mut[path$node])
#                   if (length(past_muts) > 0) past_muts[length(past_muts)] else ''
#                 },
#                 .mode = "out"
#         ),
#         FUN = identity,
#         FUN.VALUE = character(1)
#       )
#     ) %>%
#     select(-.root) %>%
#     mutate(GT = ifelse(GT == '' & !is.na(site), site, GT))
#
#   # preserve the clone ids
#   if ('GT' %in% colnames(mut_nodes)) {
#     gtree = gtree %>% activate(nodes) %>%
#       left_join(
#         mut_nodes %>% select(GT, clone),
#         by = 'GT'
#       )
#   }
#
#   return(gtree)
# }
#
#
#

#' Score the centroids tip
#' @description Score the centroids tips (can regarded as pseudo-internal node in cell-level
#' phylogeny, but is actual tips in our centroid-level phylogeny)
#' @param P matrix Genotype probability matrix
#' @param named_clusters named vector Cell assignments to each cluster
#' @return tbl_graph Phylogeny annotated with branch lengths and mutation events
score_tree_centroids_tips <- function(P, named_clusters) {

  # get q_ij for all cells
  # logQ is a (cells x sites) matrix
  logP_0 = log(1-P)
  logP_1 = log(P)

  logQ = logP_1 - logP_0


  # q_zj = sum of q_ij for all cells i belongs to cluster z
  grp  <- named_clusters[rownames(P)]

  # sum rows by cluster
  # logQ_centroids is a (clusters x sites) matrix
  logQ_centroids <- rowsum(logQ, group = grp, reorder = FALSE)

  rownames(logQ_centroids) <- unique(grp)
  colnames(logQ_centroids) <- colnames(P)

  return(logQ_centroids)
}




#' Get a tidygraph tree with simplified mutational history.
#' @description Specify either max_cost or n_cut.
#' max_cost works similarly as h and n_cut works similarly as k in stats::cutree.
#' The top-level normal diploid clone is always included.
#' @param tree phylo centroids-level phylogenetic tree
#' @param P matrix Genotype probability matrix
#' @param logQ_centroids logQ for tips (centroids) in centroids-level phylogeny
#' @param max_cost numeric Likelihood threshold to collapse internal branches
#' @param n_cut integer Number of cuts on the phylogeny to define subclones
#' @return tbl_graph Phylogeny annotated with branch lengths and mutation events
#' @export
get_gtree_centroids = function(tree, P, named_clusters, logQ_centroids, n_cut = 0, max_cost = 0) {
  # tree should be a centroids-level phylogeny
  # calculate L matrix and create tidygraph object
  tree_post = get_tree_post_centroids(tree, P, logQ_centroids)
  # print(colnames(P))
  #
  # plot(tree_post$gtree, show.tip.label = TRUE)
  # lm <- as.data.frame(tree_post$l_matrix)
  # print(lm)
  # stopifnot(identical(colnames(lm), colnames(P)))    # columns (sites) must match P
  #
  # winners <- rownames(lm)[apply(lm, 2, which.max)]
  # cat("Node winners per site (table):\n"); print(table(winners))
  # # simplify mutational history
  # nd <- tree_post$gtree %>% tidygraph::activate(nodes) %>% as.data.frame()
  # print(table(nd$last_mut))
  G_m <- get_mut_graph_centroids(tree_post$gtree)

  if (igraph::ecount(G_m) > 0) {
    G_m <- simplify_history(G_m, tree_post$l_matrix, max_cost = max_cost, n_cut = n_cut)
  }

  G_m <- label_genotype(G_m)

  # G_m = get_mut_graph_centroids(tree_post$gtree)  %>%
  #   simplify_history(tree_post$l_matrix, max_cost = max_cost, n_cut = n_cut) %>%
  #   label_genotype()

  mut_nodes = G_m %>% igraph::as_data_frame('vertices') %>%
    select(name = node, site = label, clone = clone, GT = GT)

  # update tree
  gtree = mut_to_tree(tree_post$gtree, mut_nodes)
  # attach number of cells of each centroid as a column to gtree nodes
  gtree = attach_centroid_counts(gtree, named_clusters)
  gtree = mark_tumor_lineage_centroids(gtree)

  return(gtree)

}

# named_clusters: a named vector, names = cell IDs, values = centroid labels
# centroid labels match tree tip labels used by get_gtree_centroids()

attach_centroid_counts <- function(gtree, named_clusters) {
  # create a 2-column table
  # col1 - centroids label, col2 - number of cells in that centroid
  counts <- as.data.frame(table(named_clusters), stringsAsFactors = FALSE)
  colnames(counts) <- c("name", "n_cells")

  gtree %>%
    tidygraph::activate(nodes) %>%
    dplyr::left_join(counts, by = "name")
}


#' Find maximum lilkelihood assignment of mutations on a tree
#' @param tree phylo centroids-level phylogenetic tree
#' @param P matrix Genotype probability matrix (cells x sites)
#' @param logQ_centroids logQ for tips (centroids) in centroids-level phylogeny
#' @return list Mutation
#' @keywords internal
get_tree_post_centroids = function(tree, P, logQ_centroids) {

  # tree should be a centroids-level phylogeny
  # P is still cells x sites matrix
  sites = colnames(P)
  # n = nrow(P)

  tree_stats = score_tree_centroids(tree, P, logQ_centroids, get_l_matrix = TRUE)

  # l_matrix: (tree.Ntips + tree.Nnodes) x sites
  # tree should be a centroids-level phylogeny
  l_matrix = as.data.frame(tree_stats$l_matrix)

  colnames(l_matrix) = sites
  rownames(l_matrix) = c(tree$tip.label, paste0('Node', 1:tree$Nnode))

  gtree = annotate_tree_centorids(tree, P, logQ_centroids)

  return(list('gtree' = gtree, 'l_matrix' = l_matrix))
}


#' Find maximum lilkelihood assignment of mutations on a tree
#' @param tree phylo Single-cell phylogenetic tree
#' @param P matrix Genotype probability matrix
#' @return tbl_graph A single-cell phylogeny with mutation placements
#' @examples
#' gtree_small = annotate_tree(tree_small, P_small)
#' @export
annotate_tree_centorids = function(tree, P, logQ_centroids) {

  # tree should be a centroids-level phylogeny
  sites = colnames(P)
  n = Ntip(tree)
  tree_stats = score_tree_centroids(tree, P, logQ_centroids, get_l_matrix = TRUE)

  l_matrix = as.data.frame(tree_stats$l_matrix)

  colnames(l_matrix) = sites
  rownames(l_matrix) = c(tree$tip.label, paste0('Node', 1:tree$Nnode))

  # mutation assignment on nodes
  mut_nodes = data.frame(
    site = sites,
    node_phylo = apply(l_matrix, 2, which.max),
    l = apply(l_matrix, 2, max)
  ) %>%
    mutate(
      name = ifelse(node_phylo <= n, tree$tip.label[node_phylo], paste0('Node', node_phylo - n))
    ) %>%
    group_by(name) %>%
    summarise(
      site = paste0(sort(site), collapse = ','),
      n_mut = n(),
      l = sum(l),
      .groups = 'drop'
    )

  gtree = tree %>%
    ladderize() %>%
    as_tbl_graph() %>%
    mutate(
      leaf = node_is_leaf(),
      root = node_is_root(),
      depth = bfs_dist(root = 1),
      id = row_number()
    )

  # leaf annotation for edges
  gtree = gtree %>%
    activate(edges) %>%
    select(-any_of(c('leaf'))) %>%
    left_join(
      gtree %>%
        activate(nodes) %>%
        data.frame() %>%
        select(id, leaf),
      by = c('to' = 'id')
    )

  # annotate the tree
  # we don't need to change mut_to_tree function
  gtree = mut_to_tree(gtree, mut_nodes)

  return(gtree)
}


#' Score a tree based on maximum likelihood
#' @param tree phylo object centroids-level phylogeny
#' @param P genotype probability matrix (cells x sites)
#' @param logQ_centroids logQ for tips (centroids) in centroids-level phylogeny
#' @param get_l_matrix whether to compute the whole likelihood matrix
#' @return list Likelihood scores of a tree
#' @examples
#' tree_likelihood = score_tree(tree_upgma, P_small)$l_tree
#' @export
score_tree_centroids = function(tree, P, logQ_centroids, get_l_matrix = FALSE) {

  # tree is binary from hc, so we should be safe here.
  tree = reorder(tree, order = 'postorder')

  n = Ntip(tree)
  m = ncol(P)

  # tree is binary from hc, so we should be safe here.
  logQ = matrix(nrow = tree$Nnode * 2 + 1, ncol = m)

  logP_0 = log(1-P)
  # logP_1 = log(P)

  node_order = c(tree$edge[,2], n+1)
  node_order = node_order[node_order > n] # filter for internal nodes

  logQ[1:n,] = logQ_centroids[tree$tip.label, , drop = FALSE]

  children_dict = allChildrenCPP(tree$edge)

  logQ = CgetQ(logQ, children_dict, node_order)

  if (get_l_matrix) {
    l_matrix = sweep(logQ, 2, colSums(logP_0), FUN = '+')
    l_tree = sum(apply(l_matrix, 2, max))
  } else {
    l_matrix = NULL
    l_tree = sum(apply(logQ, 2, max)) + sum(logP_0)
  }

  # l_matrix: (tree.Ntips + tree.Nnodes) x sites
  # logQ: q_v matrix for centroids-level phylogeny, (tree.Ntips + tree.Nnodes) x sites
  return(list('l_tree' = l_tree, 'logQ' = logQ, 'l_matrix' = l_matrix))

}



#' Mark the tumor lineage of a phylogeny
#' @param gtree tbl_graph Single-cell phylogeny
#' @return tbl_graph Phylogeny annotated with tumor versus normal compartment
#' @keywords internal
mark_tumor_lineage_centroids = function(gtree) {

  mut_nodes = gtree %>%
    activate(nodes) %>%
    filter(!is.na(site)) %>%
    as.data.frame() %>%
    pull(id)

  mut_burdens = lapply(
    mut_nodes,
    function(node) {
      gtree %>%
        activate(nodes) %>%
        mutate(
          mut_burden = ifelse(GT == '', 0, str_count(GT, ',') + 1)
        ) %>%
        ungroup() %>%
        mutate(seq = bfs_rank(root = node)) %>%
        data.frame %>%
        filter(leaf & seq > 0) %>%
        mutate(w = n_cells) %>%
        mutate(total = w * mut_burden) %>%
        pull(total) %>%
        # pull(mut_burden) %>%
        sum
    }
  )

  tumor_root = mut_nodes[which.max(mut_burdens)]

  gtree = gtree %>%
    activate(nodes) %>%
    mutate(
      seq = bfs_rank(root = tumor_root),
      compartment = ifelse(seq > 0, 'tumor', 'normal'),
      is_tumor_root = tumor_root == id
    )

  compartment_dict = gtree %>%
    activate(nodes) %>%
    as.data.frame() %>%
    {setNames(.$compartment, .$id)}

  gtree = gtree %>%
    activate(edges) %>%
    mutate(compartment = compartment_dict[to])

  return(gtree)

}

#' Convert a centroids-level phylogeny with mutation placements into a mutation graph
#'
#' @param gtree tbl_graph The phylogeny
#' @return igraph Mutation graph
#' @examples
#' mut_graph = get_mut_graph_centroids(gtree_small)
#' @export
get_mut_graph_centroids = function(gtree) {

  mut_nodes = gtree %>%
    activate(nodes) %>%
    as.data.frame() %>%
    dplyr::filter(!is.na(site)) %>%
    dplyr::distinct(name, site)

  # if last_mut only defines one group, return a single-vertex mutation graph
  node_df = gtree %>%
    activate(nodes) %>%
    as.data.frame()

  # check last_mut
  lm = node_df$last_mut
  if (is.null(lm)) {
    lm = rep(NA_character_, nrow(node_df))
  }

  u_last = unique(lm[!is.na(lm)])

  if (length(u_last) <= 1) {

    # build a 1-vertex igraph carrying the available label ('' or the only last_mut)
    lab = if (length(u_last) == 1) u_last[[1]] else ''

    G = igraph::make_empty_graph(n = 1, directed = TRUE)
    igraph::V(G)$label = lab
    igraph::V(G)$id    = 1

    # try to map back to a node name if the label matches a site in mut_nodes
    node_map = mut_nodes %>% dplyr::rename(node = name)
    m = node_map$node[match(lab, node_map$site)]
    igraph::V(G)$node = ifelse(is.na(m), NA_character_, m)

    return(G)
  }


  # if last_mut only defines 2 or more group, return a mutation graph as original Numbat
  G = gtree %>%
    activate(nodes) %>%
    dplyr::arrange(last_mut) %>%
    convert(to_contracted, last_mut) %>%
    dplyr::mutate(label = last_mut, id = 1:n()) %>%
    as.igraph

  G = label_edges(G)

  igraph::V(G)$node = G %>%
    igraph::as_data_frame('vertices') %>%
    dplyr::left_join(
      mut_nodes %>% dplyr::rename(node = name),
      by = c('label' = 'site')
    ) %>%
    dplyr::pull(node)

  G = G %>% transfer_links()

  return(G)
}

