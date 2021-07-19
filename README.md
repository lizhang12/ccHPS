## ccHPS
a Hypoxia-related Prognostic Signature for cervical cancer (ccHPS), which was developed and validated 
using tumor samples in TCGA-CESC and CGCI-HTMCP-CC cohorts based on LASSO Cox regression analysis.

## ccHPS.R
calculate the ccHPS score for each input samples using gene expression matrix.
# Usage:
```
Rscript ccHPS.R expression_input_file.txt
```
`input` The gene expression matrix file is fomated as genes in each column and samples in each row. 
The input genes should include all the nine genes involved in ccHPS model.
An example format of `expression_input_file.txt`:

![Uploading Screenshot 2021-07-19 181327.png…]()

`output` an output file will be created with the ccHPS score of each sample.


