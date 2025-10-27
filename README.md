# tbcluster-basic-snp


# Libraries

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
vcf_file <- "D:/Bioinformatics/Thailand_TB_Data/L2.2.M3_HAP1-2/snpplet/results/joint_called/filtered_vcf/joint_filtered.vcf.gz 
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

```
cutoff <- 12
g <- graph.adjacency(snp_mat_for_graph <= cutoff, mode="undirected", diag=FALSE)
clusters <- components(g)$membership

cluster_table <- data.frame(Sample=names(clusters), Cluster=clusters)
table(clusters)
head(cluster_table)

```

# 7️⃣ Visualization using ggraph

```
n_clusters <- max(clusters)
cluster_colors <- brewer.pal(min(n_clusters,12), "Set3")
V(g)$cluster <- clusters
V(g)$color <- cluster_colors[clusters]

set.seed(123)
ggraph(g, layout="fr") +
  geom_edge_link(color="gray50") +
  geom_node_point(aes(color=factor(cluster)), size=5) +
  theme_void() +
  scale_color_brewer(palette="Set3") +
  ggtitle(paste("TB clusters with", cutoff, "SNP cutoff"))

```
