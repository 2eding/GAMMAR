# GRAMMAR
Generalized Renown Analysis of Molecular variance for Mixed model Analysis R package

# usage
library(data.table)

X = as.matrix(fread("X_rightdim.txt"))
Y = as.matrix(fread("Y_rightdim.txt"))

K = Kinship(t(X))
VC = varComp(K, Y, X)
run = run_grammar(VC[3], VC[4], max_itr = 4, num.parallel = 4)