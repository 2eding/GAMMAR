# GRAMMAR
Generalized Renown Analysis of Molecular variance for Mixed model Analysis R package

# usage
library(data.table)<br>
library(devtools)<br>
install_github("2eding/GRAMMAR")<br>
library(GRAMMAR)<br><br>

X = as.matrix(fread("X_rightdim.txt"))<br>
Y = as.matrix(fread("Y_rightdim.txt"))<br><br>

K = Kinship(X)<br>
VC = varComp(Y, K, num.parallel)<br>
run = run_grammar(K, Y, X, VC, max_itr = 4, num.parallel = 4, outPath = "./")
