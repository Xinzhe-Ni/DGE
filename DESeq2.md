# DGE with DESeq2

在基因差异表达分析中，最常用的R包为DESeq2和edgeR。这里将用DESeq2来演示如何进行基因差异表达分析，用到的数据集链接在这里。这里对野生型（WT）在光照前和光照后的基因表达差异进行演示。

* **数据说明**

  > 该数据集是拟南芥野生型（WT）与uvr8突变型（*uvr8*）在光照前后的基因表达矩阵。第1列为基因名称，第2～4列为WT光照前的表达值，第5～7列为WT光照后的表达值，第8～10列为*uvr8*光照前的表达值，第11～13列为*uvr8*光照后的表达值。这里演示只用到前7列的数据。

## 1\) Input data

``` R
library(DESeq2)

raw_count <- read.table("count_exon.txt", sep='\t', header = T)

#提取野生型
wt_raw_count <- raw_count[c("gene_id", "CD1_1", "CD1_2", "CD1_3", "CD0_1", "CD0_2", "CD0_3")]
row.names(wt_raw_count) <- wt_raw_count[, 1]
wt_raw_count <- wt_raw_count[, -1]

#过滤
countData_wt <- wt_raw_count[rowSums(wt_raw_count) > 100, ]
```

第一步：读取数据，得到野生型的基因表达矩阵。

> Tips: 输入矩阵需为原始矩阵，不能经过标准化等处理

## 2\) Provide sample information

``` R
condition_merge <- factor(c(rep("Control", 3), rep("Treat", 3)))
colData <- data.frame(row.names = colnames(countData_wt), condition_merge)
```

第二步：提供分组信息，光照前与光照后记作Control和Treat组，记录在信息矩阵中。

## 3\) DGE analysis

``` R
#获取dds矩阵
dds <- DESeqDataSetFromMatrix(countData_wt, colData, design = ~condition_merge)

#标准化
dds2 <- DESeq(dds)

#获取结果
res <- results(dds2)
```

DESeq包含三步：estimation of size factor (estimateSizeFactors) , estimation of dispersion (estimateDispersons) , Negative Binomial GLM fitting and Wald statistics (nbinomWaldTest)。

## 4\) Filter

``` R
res <- res[order(res$padj), ]

#过滤标准
diff_gene_deseq_wt <- subset(res, padj < 0.05 & abs(log2FoldChange > 1))

#提取差异基因名称
wt_gene_names <- row.names(diff_gene_deseq_wt)
write.csv(wt_gene_names,"wt_gene_list.txt", sep='\t', row.names = F, quote = F)
```

提取基因名称之后，可进行下一步如GO、KEGG的分析。

## 5\) Homework

按照教程中的流程，使用DESeq2找出uvr8突变型（*uvr8*）在光照前后的差异基因，输出一个txt文件。

