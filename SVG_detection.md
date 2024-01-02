SVG_detection
================
HSY
2024-01-01

> 运行环境中的版本信息如下：R (4.1.0), Seurat (4.0.0), SPARK (1.1.1)

## SPARK

``` r
library(SPARK)
library(Seurat)

oneseuv = readRDS("data/ST_Seurat.rds")
rawcount = oneseuv@assays$Spatial@counts

## annotation information for each sample
info = oneseuv@images$slice1@coordinates[,c("imagecol","imagerow")]
colnames(info) = c("x","y")

## Create a SPARK object for analysis
## filter genes and cells/spots
spark <- CreateSPARKObject(
  counts=rawcount, 
  location=info,
  percentage = 0.05, #取值范围0-1，某一个gene表达的spot数量必须>=（percentage * 总spot数量）
  min_total_counts = 10 #每个有效的spot中，必须大于这个umi数
)

## total counts for each cell/spot
spark@lib_size <- apply(spark@counts, 2, sum)

## 为了演示，仅选取少量基因，其中OLFM4明确是svg（后来的分析结果表明的）
example_gene = head(rownames(spark@counts))
example_gene = c(example_gene,"OLFM4")
spark@counts   <- spark@counts[example_gene,]

## Fit the statistical model under the null hypothesis
## Estimating Parameter Under Null
spark <- spark.vc(
  spark, 
  covariates = NULL, 
  lib_size = spark@lib_size, 
  num_core = 1,
  verbose = T)

## Test the spatially expressed pattern genes
## Calculating pval
## 这一步较慢
spark <- spark.test(
  spark, 
  check_positive = T, 
  verbose = T)

res = spark@res_mtest[,c("combined_pvalue","adjusted_pvalue")]
res

saveRDS(spark,file = "data/spark.rds")
```

## SPARK-X

适用于spot数量较多的情况

``` r
library(SPARK)
library(Seurat)
```

    ## Registered S3 method overwritten by 'spatstat':
    ##   method     from
    ##   print.boxx cli

    ## Attaching SeuratObject

``` r
oneseuv = readRDS("data/ST_Seurat.rds")
rawcount = oneseuv@assays$Spatial@counts
info = oneseuv@images$slice1@coordinates[,c("imagecol","imagerow")]
colnames(info) = c("x","y")

sparkX <- sparkx(rawcount,info,numCores=1,option="mixture")
```

    ## ## ===== SPARK-X INPUT INFORMATION ==== 
    ## ## number of total samples: 3138 
    ## ## number of total genes: 23233 
    ## ## Running with single core, may take some time 
    ## ## Testing With Projection Kernel
    ## ## Testing With Gaussian Kernel 1
    ## ## Testing With Gaussian Kernel 2
    ## ## Testing With Gaussian Kernel 3
    ## ## Testing With Gaussian Kernel 4
    ## ## Testing With Gaussian Kernel 5
    ## ## Testing With Cosine Kernel 1
    ## ## Testing With Cosine Kernel 2
    ## ## Testing With Cosine Kernel 3
    ## ## Testing With Cosine Kernel 4
    ## ## Testing With Cosine Kernel 5

``` r
res.sx <- sparkX$res_mtest
head(res.sx)
```

    ##            combinedPval adjustedPval
    ## AL627309.1 8.470885e-01  1.000000000
    ## AL627309.5 9.826298e-01  1.000000000
    ## AP006222.2 2.283793e-01  1.000000000
    ## AL669831.2 8.908458e-01  1.000000000
    ## LINC01409  9.155142e-02  1.000000000
    ## LINC01128  6.890367e-06  0.000234923
