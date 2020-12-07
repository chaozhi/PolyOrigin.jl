# Installing and using the PedigreeSimR to simulate the cross and GBS data
# devtools::install_github("rramadeu/PedigreeSimR",
#                          auth_token = "47b8aa6fe594e5c42efc7f4c0b951a2efc0b6663",
#                          force=TRUE)
library(PedigreeSimR)
packageVersion("PedigreeSimR")

library(rstudioapi)    
workdir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
getwd()

simmaphaplo <- function(nsnp,np,ploidy=4,chrlen=100,seed=1234){
  haplo = fake_haplo(n=np*ploidy,m=nsnp,seed=seed)
  cm0=rexp(nsnp-1,rate=1)
  cm0=cm0/(sum(cm0)/chrlen)
  cm0=round(cm0,digits=2)
  cm=c(0,cumsum(cm0))
  list(map=cm,hap=haplo)
}

ploidy=4
nsnp=120
np=3
sim=simmaphaplo(nsnp,np,ploidy=ploidy)
pedigree = diallel_pedigree(parents=np,selfs=1,popsize=40)
pedigreesimR(sim$map,sim$hap,
             ploidy=ploidy,
             sampleHap = FALSE,filename="4x_diallel_",
             pedigree,workingfolder = getwd(),
             epsilon = c(0.01,0.01),
             missingFreq=c(0.1,0.1),
             # GBS=TRUE,GBSavgdepth = 20
)

# delete all pedsim*.* files. 
sapply(list.files(pattern = "pedsim*"), unlink)
unlink("*truevalue.csv")
unlink("*geno.csv")
