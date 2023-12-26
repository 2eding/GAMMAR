# GAMMAR
Generalized Renown Analysis of Molecular variance for Mixed model Analysis R package

# Usage

### Step 1. Install and load the package
library(devtools)<br>
install_github("2eding/GAMMAR")<br>
library(GAMMAR)<br><br>

### Step 2. Input data
X <- as.matrix(data.table::fread("X_rightdim.txt"))<br>
Y <- as.matrix(data.table::fread("Y_rightdim.txt"))<br>
* data.table::fread <= This package is useful for input large data<br>
* X_rightdim.txt is snp data
* Y_rightdim.txt is gene expression data
<br><br>
### Step 3. Calculate the kinship
K <- Kinship(X)<br><br>

### Step 4. Calculate the variance components
VC <- varComp(K, Y, X)<br><br>

### Step 5. Run GAMMAR
run <- run_grammar(K, Y, X, VC, max_itr = 4, num.parallel = 4, outPath = "./test/", outName = "result.txt")

# Publication
http://doi.org/10.5391/IJFIS.2020.20.3.1
