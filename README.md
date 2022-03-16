## ccHPS
a Hypoxia-related Prognostic Signature for cervical cancer (ccHPS), which was developed and validated 
using tumor samples in TCGA-CESC and CGCI-HTMCP-CC cohorts based on LASSO Cox regression analysis. The data used in the developing process could be found in this link https://figshare.com/articles/dataset/Data_set_repositories/14999487.

## ccHPS.R
calculate the ccHPS score for each input samples using gene expression matrix.
# Usage:
```
Rscript ccHPS.R expression_input_file.txt
```
`input` The gene expression matrix file is fomated as genes in each column and samples in each row. 
The input genes should include all the nine genes involved in ccHPS model.
An example format of `expression_input_file.txt`:

![Screenshot 2021-07-19 181327](https://user-images.githubusercontent.com/61243260/126146243-65cd1150-f991-45b3-91fb-76fda262588a.png)


`output` an output file will be created with the ccHPS score of each sample.


