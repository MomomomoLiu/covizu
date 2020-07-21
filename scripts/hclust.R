require(igraph)
require(jsonlite)
require(Rtsne)
require(ape)

# open TN93 distance matrix
cat("loading TN93 distance matrix\n")
tn93 <- read.csv('data/variants.tn93.txt', skip=1, header=F)
tn93 <- as.matrix(tn93)
stopifnot(nrow(tn93) == ncol(tn93))


# apply hierarchical clustering to distance matrix
cat("hierarchical clustering\n")

# direct clustering
# hc <- hclust(as.dist(tn93), method='complete')
# clusters <- cutree(hc, h=0.002)
# hist(log10(table(clusters)), breaks=20, col='grey', border='white')

# cluster t-stochastic neighbor embedding
set.seed(1)
res <- Rtsne(tn93, is_distance=TRUE, verbose=TRUE, num_threads=2)
hc <- hclust(dist(res$Y))
clusters <- cutree(hc, h=10)

# use regex to extract headers from FASTA
headers <- rep(NA, times=nrow(tn93))
con <- file('data/variants.fa', open='r')
i <- 1
L <- 0
while (length(line <- readLines(con, n=1, warn=FALSE)) > 0) {
  if (grepl("^>", line)) {
    headers[i] <- gsub("^>(.+)", "\\1", line)
    i <- i + 1
  }
  else if (L == 0) {
    L <- nchar(line)
  }
}
close(con)
stopifnot(length(headers) == nrow(tn93))
names(tn93) <- headers
accns <- gsub("^.+(EPI_[A-Z]+_[0-9]+).+$", "\\1", headers)

# extract sample collection dates from headers
dates <- as.Date(sapply(headers, function(x) {
  strsplit(x, "\\|")[[1]][3]
  }))

# identify the earliest sample (Wuhan, IPBCAMS-WH-01)
#root <- which.min(dates)
root <- which(grepl("EPI_ISL_402132", accns))


tab <- table(clusters)
mean.pdist <- c()
mean.rdist <- c()
for (k in 1:length(tab)) {
  clust <- as.integer(which(clusters == k))

  # compute mean pairwise distance between members of each cluster  
  pdists <- as.matrix(tn93[clust, clust])
  mdist <- mean(pdists[upper.tri(pdists)])
  mean.pdist <- c(mean.pdist, mdist)

  # compute mean distance between members and the root  
  rdists <- as.matrix(tn93[root, clust])
  mrdist <- mean(rdists)
  mean.rdist <- c(mean.rdist, mrdist)
}

# open CSV with SARS-COV-2 genome variant information
variants <- read.csv('data/variants.csv')
stopifnot(all(is.element(variants$cluster, headers)))

#' @param node: str, label of current node variant
#' @param parent: str, label of current node's parental variant
#' @param el: str, edge list from minimum spanning tree
#' @return linearized vector of parent->child pairs
traverse <- function(node, parent, el, edges=c()) {
  if (!is.na(parent)) {
    edges <- c(edges, parent, node)  
  }
  # get adjacent to current node
  temp <- el[apply(el, 1, function(e) is.element(node, e)), ]
  temp <- unique(as.vector(temp))
  children <- temp[!is.element(temp, c(node, parent))]
  for (child in children) {
    edges <- traverse(child, node, el, edges)
  }
  return(edges)
}


count.tips <- function(edge.mx) {
  # convert matrix to data frame
  edge.df <- as.data.frame(edge.mx, stringsAsFactors = FALSE)
  names(edge.df) <- c('parent', 'child')
  edge.df$n.tips <- 0
  edge.df$coldate <- sapply(edge.df$child, function(x) {
    as.Date(strsplit(x, "\\|")[[1]][3])
  })
  
  # retrieve row indices of tips
  idx <- !is.element(edge.df$child, edge.df$parent)
  edge.df$n.tips[idx] <- 1
  children <- edge.df$child[idx]  # tip labels
  
  while (length(children) > 0) {
    # retrieve parents
    idx <- match(children, edge.df$child)
    parents <- edge.df$parent[idx]
    
    # find next parents
    idx2 <- unique(match(parents, edge.df$child))
    
    temp <- sapply(
      split(idx, parents), function(x) sum(edge.df$n.tips[x])
    )
    
    if (any(is.na(idx2))) {
      drop <- which(is.na(idx2))
      idx2 <- idx2[-drop]
      temp <- temp[-drop]
    } 
    
    # propagate tip counts
    edge.df$n.tips[idx2] <- temp
    children <- edge.df$child[idx2]
  }
  
  edge.df
}


reorder.edges <- function(edge.df, tn93) {
  # augment data frame with TN93 distances
  edge.df$dist <- sapply(1:nrow(edge.df), function(i) {
    row <- match(edge.df$parent[i], colnames(tn93))
    col <- match(edge.df$child[i], colnames(tn93))
    tn93[row, col]
  })
  
  idx <- match(edge.df$parent, unique(edge.df$parent))
  idx2 <- as.vector(unlist(lapply(split(1:nrow(edge.df), idx), function(i) {
    rows <- edge.df[i, ]
    i[order(rows$n.tips, rows$dist, rows$coldate)]
  })))
  edge.df[idx2, ]
}


result <- list()
for (i in 1:max(clusters)) {
  cat('.')  # crude progress monitoring
  
  # extract cluster indices to map to headers vector
  idx <- as.integer(which(clusters==i))

  pdist <- mean.pdist[[i]]  # cluster mean pairwise distance
  rdist <- mean.rdist[[i]]  # cluster mean root distance
  
  if (length(idx)==1) {
    list(nodes=headers[idx], edges=NA, pdist=pdist, rdist=rdist)
  } 
  else {
    # find earliest variant in cluster that is closest to root
    min.dist <- min(tn93[root, idx])
    candidates <- headers[idx[which(tn93[root, idx] == min.dist)]]
    subroot <- candidates[which.min(sapply(candidates, function(x) { 
      as.Date(strsplit(x, "\\|")[[1]][3])
      }))]
    
    # generate minimum spanning tree
    mx <- as.matrix(tn93[idx, idx])
    colnames(mx) <- headers[idx]
    g <- graph.adjacency(mx, weighted=TRUE)
    g.mst <- igraph::mst(g, algorithm='prim')
    
    # traverse MST and export node and edge lists
    el <- get.edgelist(g.mst)
    edges <- traverse(subroot, NA, el)

    edge.mx <- matrix(edges, ncol=2, byrow=TRUE)
    edge.df <- count.tips(edge.mx)
    edge.df <- reorder.edges(edge.df, tn93)
    
    # store variant data
    nodes <- list()
    for (node in unique(edges)) {
      accn <- strsplit(node, "\\|")[[1]][2]
      temp <- variants[variants$cluster==node, ]
      temp$label1 <- sapply(as.character(temp$label), function(x) {
        strsplit(x, "\\|")[[1]][1]
      })
      temp$accession <- sapply(as.character(temp$label), function(x) {
        strsplit(x, "\\|")[[1]][2]
      })
      nodes[[accn]] <- temp[c('label1', 'accession', 'country', 'coldate')]
    }

    # shorten edge list to accession numbers only
    edges <- cbind(
      gsub("^.+(EPI_[A-Z]+_[0-9]+).+$", "\\1", edge.df$parent),
      gsub("^.+(EPI_[A-Z]+_[0-9]+).+$", "\\1", edge.df$child),
      round(edge.df$dist * L, 2)
    )
    
    result[[length(result)+1]] <- list(
      pdist=pdist, rdist=rdist, nodes=nodes, edges=edges
      )
  }
}
cat ('\nwriting JSON file\n')
write(toJSON(result, pretty=TRUE), file="data/clusters.json")

