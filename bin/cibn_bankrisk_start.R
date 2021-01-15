mainDir <- "D:/Research/R Packages/cibn"


source(file.path(mainDir,"R/cibn.r"))

# create the BN
bankrisk_bn <- new.cibn(path=file.path(mainDir,"bin/bankrisk_bn.txt"))

# plot the DAG
plot(bankrisk_bn, attrs=list(edge=list(arrowsize=0.5)))
#plot(bankrisk_bn, full=T, attrs=list(edge=list(arrowsize=0.5)))

# evidence propagation
getStates(bankrisk_bn)
query.cibn(bankrisk_bn, target="Risk", evidence=list(Age="31_50",Portf="mixed"))
query.cibn(bankrisk_bn, target="Risk", evidence=list(Age="31_50",Portf="money_market"))
query.cibn(bankrisk_bn, target="Risk", evidence=list(Age="31_50",Portf="stock_market"))
