# gsea_visualisation

R object of example gsea results is here:
http://ftp1.babraham.ac.uk/sites/snappy-snail-8241/ 

Link to gene set course exercises:
https://www.bioinformatics.babraham.ac.uk/training/Gene_Set_Analysis/R_Gene_Set_Analysis_Exericse.pdf 

The code that was run to create the gsea_results object (it takes quite a long time to run):

(`ranks` is the ranked_genes.rds object in this repo)

```
gsea_results  <- gseGO(
  geneList     = ranks,
  ont          = "ALL",
  OrgDb        = org.Hs.eg.db,
  keyType      = "SYMBOL",
  pvalueCutoff = Inf,
  minGSSize    = 10,
  maxGSSize    = 200
)
```

The second gsea_results object on the sftp site was generated using this code:

```
library(DOSE)
data(geneList)

ranks <- geneList

gsea_results  <- gseGO(
  geneList     = ranks,
  ont          = "ALL",
  OrgDb        = org.Hs.eg.db,
  keyType      = "ENTREZID",
  pvalueCutoff = Inf,
  minGSSize    = 10,
  maxGSSize    = 200
)
```
