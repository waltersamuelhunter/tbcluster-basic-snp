
# Packages
```
install.packages(c(
  "vcfR",
  "igraph",
  "ggraph",
  "ggplot2",
  "RColorBrewer",
  "tidygraph",
  "ggmap",
  "scales",
  "readxl",
  "devtools"
  ))

devtools::install_github("JamesStimson/transcluster", build_vignettes = TRUE)

```
Load packages 
```
library(vcfR)
library(igraph)
library(ggraph)
library(ggplot2)
library(RColorBrewer)
library(tidygraph)
library(ggmap)
library(scales)
library(readxl)

```

# 1️⃣ Read VCF

Read the joint vcf file - specify file location!

```
vcf_file <- "D:/specify/your/path/.vcf.gz"
vcf_r <- read.vcfR(vcf_file)
head(vcr_r)
```

# 2️⃣ Filter variants and samples with too much missing data

```
geno_mat <- extract.gt(vcf_r, element="GT", return.alleles = TRUE)  # rows = variants, cols = samples
```

# Filter variants: keep if <= 20% missing (try different values) 

```
variant_na_frac <- apply(geno_mat, 1, function(x) mean(is.na(x) | x == "." | x == "./."))
ggplot(data=data.frame(x=variant_na_frac))+
  geom_histogram(aes(x=x), bins = 100)+
  labs(x="Fraction",y="Variant Count")
```
```
keep_variants <- variant_na_frac <= 0.2
length(variant_na_frac)
geno_mat <- geno_mat[keep_variants, ]
```

# Filter samples: keep if <= 20% missing

```
sample_na_frac <- apply(geno_mat, 2, function(x) mean(is.na(x) | x == "." | x == "./."))
ggplot(data=data.frame(x=sample_na_frac))+
  geom_histogram(aes(x=x), bins = 150)+
  labs(x="Fraction", y = "Sample Count")
```
```
keep_samples <- sample_na_frac <= 0.25
geno_mat <- geno_mat[, keep_samples]
ggplot(data=data.frame(x=sample_na_frac))+
  geom_histogram(aes(x=x), bins = 100)
geno_num = geno_mat
#geno_num <- apply(geno_mat, 2, function(x) as.numeric(x))
```


# 3️⃣ Pairwise SNP distance ignoring NAs

```
pairwise_snp <- function(x, y) {
  valid <- !is.na(x) & !is.na(y)
  if(sum(valid) == 0) return(NA)
  sum(x[valid] != y[valid])
}

samples <- colnames(geno_num)
n_sam <- length(samples)
snp_mat <- matrix(NA, nrow=n_sam, ncol=n_sam, dimnames=list(samples, samples))

for(i in 1:n_sam) {
  for(j in i:n_sam) {
    d <- pairwise_snp(geno_num[,i], geno_num[,j])
    snp_mat[i,j] <- d
    snp_mat[j,i] <- d
  }
}

```

# 4️⃣ Prepare matrix for clustering

```
snp_mat_for_graph <- snp_mat
snp_mat_for_graph[is.na(snp_mat_for_graph)] <- 10000  # disconnected
```

# 5️⃣ 12-SNP cutoff clustering using igraph 

You can adjust the cutoff value, to change clustering parameters.

```
cutoff <- 12
g <- graph.adjacency(snp_mat_for_graph <= cutoff, mode="undirected", diag=FALSE)
clusters <- components(g)$membership

cluster_table <- data.frame(Sample=names(clusters), Cluster=clusters)
table(clusters)
head(cluster_table)

```

# 6️⃣ Visualization using ggraph

#You may play with different colour palletes by changing the value of scale_color_brewer(palette="Set3") inside ggraph 

```
n_clusters <- max(clusters)
cluster_colors <- brewer.pal(min(n_clusters,12), "Set3")
V(g)$cluster <- clusters
V(g)$color <- cluster_colors[clusters]

set.seed(123)

# Suppose you already have:
# g <- your graph
# clusters <- your clustering vector
# cutoff <- your cutoff value

# Assign clusters
V(g)$cluster <- clusters

# Identify nodes that have edges (degree > 0)
has_edges <- degree(g) > 0

# Identify which clusters actually have edges
clusters_with_edges <- unique(V(g)$cluster[has_edges])

# Assign colors only to those clusters
n_clusters <- length(clusters_with_edges)
cluster_colors <- brewer.pal(min(n_clusters, 12), "Set3")

# Map cluster colors
V(g)$color <- "grey80"  # default color for isolated nodes
color_map <- setNames(cluster_colors, clusters_with_edges)
V(g)$color[has_edges] <- color_map[as.character(V(g)$cluster[has_edges])]
E(g)$snp <- snp_mat_for_graph[as.matrix(get.edges(g, 1:ecount(g)))]

# Basic Plot
set.seed(123)

ggraph(g, layout = "fr") +
  geom_edge_link(color = "blue4", alpha = 0.8) +
  geom_node_point(aes(color = factor(cluster)), size = 4, show.legend = TRUE) +
  scale_color_manual(values = color_map,
                     name = "Clusters with edges") +
  ggtitle(paste("TB clusters with", cutoff, "SNP cutoff")) +
  theme_void()


# Plot with variable edge width, depending on SNP distance
ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = snp), color = "pink", alpha = 0.6) +
  geom_node_point(aes(color = color), size = 4) +
  scale_color_identity() + 
  scale_edge_width_continuous(breaks = seq(min(E(g)$snp), max(E(g)$snp), by = 1)) +
  ggtitle(paste("All TB samples (SNP cutoff:", cutoff, ")")) +
  theme_void()

meta <- metadata_l22m3
meta <- meta[match(V(g)$name, meta$sample_name), ]  # align rows

# Add metadata as vertex attributes
V(g)$cur_province <- meta$cur_prov

# Create a color palette for provinces
provinces <- unique(V(g)$cur_province)
province_colors <- hue_pal()(length(provinces))
province_map <- setNames(province_colors, provinces)

# Assign colors to vertices
V(g)$color <- province_map[V(g)$cur_province]

# Plot graph with province coloring
set.seed(123)
ggraph(g, layout = "fr") +
  geom_edge_link(aes(width = snp), color = "gray50", alpha = 0.6) +
  geom_node_point(aes(color = cur_province), size = 4) +
  scale_color_manual(values = province_map, name = "Province") +
  ggtitle(paste("TB clusters with metadata (SNP cutoff:", cutoff, ")")) +
  theme_void()

```

