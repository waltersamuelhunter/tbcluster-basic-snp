# tbcluster-basic-snp

The tbcluster-basic-snp is an R script for clustering Mtb genomic data. The input file is joint vcf (multiple samples) file. 

# Software 

The script requires minimal effort -  user just have to install the required packages, specify the path for joint vcf file and select the cutoff value for clustering. 

# How it works 

1. Taken all samples into account, variable sites are extracted. 
2. The samples are compared in a pairs for each analysed variant. For example, considering one variant site (v1) and two samples (x1 and x2), if sample x1 have different genotype in the v1 than sample x2 then the genetic distance between these two samples will be increased by 1. This is repeated for all samples and variant sites and the result is stored in snp distance matrix. 
3. A graph is computed containing nodes (samples) and edges (same clusters). Edges are only drawn between clusters defined by the cutoff value (like 5 or 12 SNPs)
4. The graph is visualized.


# Integration with other pipelines

For smooth analysis of multiple samples, tbcluster-basic-snp can be used in combination with nextflow pipelines. 
The joint vcf output from snpplet (https://github.com/CENMIG/snpplet) or other pipelines can be used as input from  tbcluster-basic-snp.

# Limitations 

Clustering TB cases based on SNP pairwise distance alone is known to be limited method, especially in high burden countries, due to various factors. Nevertheless, this method alone is useful for getting basic understanding of transmission in a geographical location and formulating hypothesis for public health interventions.
Also, the SNP approach can be extended to include various metadata like drug resistance, lineage, occupation etc. The script will be updated to include this options soon. 

# How to interpret

Nodes - samples 
Edges are drawn between samples (nodes), only if the snp distance between samples is less or equal to 12 SNPs. 
The connected samples can be assummed to be part of same transmission cluster.

# Packages

Install from CRAN
```
install.packages(c(
  "vcfR",
  "igraph",
  "ggraph",
  "ggplot2",
  "RColorBrewer"
))

```
Load packages 
```
library(vcfR)
library(igraph)
library(ggraph)
library(ggplot2)
library(RColorBrewer)
```

# 1️⃣ Read VCF

Read the joint vcf file - specify file location!

```
vcf_file <- "D:/specify/your/path/.vcf.gz"
vcf_r <- read.vcfR(vcf_file)

```

# 2️⃣ Filter variants and samples with too much missing data

```
geno_mat <- extract.gt(vcf_r, element="GT")  # rows = variants, cols = samples
```

# Filter variants: keep if <= 20% missing

```
variant_na_frac <- apply(geno_mat, 1, function(x) mean(is.na(x) | x == "." | x == "./."))
keep_variants <- variant_na_frac <= 0.2
geno_mat <- geno_mat[keep_variants, ]
```

# Filter samples: keep if <= 20% missing

```
sample_na_frac <- apply(geno_mat, 2, function(x) mean(is.na(x) | x == "." | x == "./."))
keep_samples <- sample_na_frac <= 0.2
geno_mat <- geno_mat[, keep_samples]
```

# 3️⃣ Convert genotypes to integers 

```
geno_num <- apply(geno_mat, 2, function(x) as.numeric(x))
```

# 4️⃣ Pairwise SNP distance ignoring NAs

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

# 5️⃣ Prepare matrix for clustering

```
snp_mat_for_graph <- snp_mat
snp_mat_for_graph[is.na(snp_mat_for_graph)] <- 10000  # disconnected
```

# 6️⃣ 12-SNP cutoff clustering using igraph 

You can adjust the cutoff value, to change clustering parameters.

```
cutoff <- 12
g <- graph.adjacency(snp_mat_for_graph <= cutoff, mode="undirected", diag=FALSE)
clusters <- components(g)$membership

cluster_table <- data.frame(Sample=names(clusters), Cluster=clusters)
table(clusters)
head(cluster_table)

```

# 7️⃣ Visualization using ggraph

You may play with different colour palletes by changing the value of scale_color_brewer(palette="Set3") inside ggraph 


```
n_clusters <- max(clusters)
cluster_colors <- brewer.pal(min(n_clusters,12), "Set3")
V(g)$cluster <- clusters
V(g)$color <- cluster_colors[clusters]

set.seed(123)
library(tidygraph)
library(igraph)
library(ggraph)
library(RColorBrewer)

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

# Plot
set.seed(123)
ggraph(g, layout = "fr") +
  geom_edge_link(color = "gray60", alpha = 0.8) +
  geom_node_point(aes(color = factor(cluster)), size = 4, show.legend = TRUE) +
  scale_color_manual(values = color_map,
                     name = "Clusters with edges") +
  ggtitle(paste("TB clusters with", cutoff, "SNP cutoff")) +
  theme_void()

```
