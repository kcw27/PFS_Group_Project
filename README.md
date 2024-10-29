# Programming For Scientists (02-601) Group Project
Course taught by Professor Phillip Compeau at Carnegie Mellon University, Fall 2024  
<a href="https://docs.google.com/document/d/1JmdNldR_mWnpVBF9NvQcTxIbkkOmzIdUMIrUF6iC4wQ/edit?usp=sharing">Google Doc</a>

## Project description
In Go, we will reimplement the experiment described in Additional File 4 of <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-11-497#Sec18">Tesson _et al._ 2010</a>. The purpose of this experiment was to demonstrate that DiffCoEx is more sensitive to differential coexpression modules than other algorithms are.

### Steps
1. Process gene expression matrices in Go, using Additional File 1 of the same paper as reference. Different algorithms may require a slightly different format for the data. Also generate toy data for testing purposes. Data may be found <a href="https://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS2901">here (uploaded 2007)</a> or <a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5923">here (updated 2017)</a>. 
2. Run a pipeline from Go to R in order to run DiffCoEx, coXpress, and potentially other algorithms on the data. We will use existing R packages for this step.
3. In an R Shiny web app, visualize the differences in output between the algorithms being compared.
It may be helpful to write executable files to automate the pipeline (e.g. automatically feed the output from our Go pre-processing application into R for analysis).

### TO DO:
1. Delegate tasks.
2. Follow the steps outlined in Additional File 1 to familiarize ourselves with the data processing expected for inputs to DiffCoEx and other algorithms.

## Citations
* Stemmer K, Ellinger-Ziegelbauer H, Ahr HJ, Dietrich DR. Carcinogen-specific gene expression profiles in short-term treated Eker and wild-type rats indicative of pathways involved in renal tumorigenesis. _Cancer Res_ 2007 May 1;67(9):4052-68. PMID: 17483316
* Tesson, B.M., Breitling, R. & Jansen, R.C. DiffCoEx: a simple and sensitive method to find differentially coexpressed gene modules. _BMC Bioinformatics_ 11, 497 (2010). https://doi.org/10.1186/1471-2105-11-497
* Watson, M. CoXpress: differential co-expression in gene expression data. _BMC Bioinformatics_ 7, 509 (2006). https://doi.org/10.1186/1471-2105-7-509
