# Installing and using the PedigreeSimR to simulate the cross and GBS data
# devtools::install_github("rramadeu/PedigreeSimR")
library(PedigreeSimR)
packageVersion("PedigreeSimR")

library(rstudioapi)    
workdir=dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(workdir)
getwd()

# set the founder haplotypes
df = read.table("TableS1_RussetFounderHaplo.csv",sep=",",header = TRUE,skip=1)
russet.map = df$cM
russet.hap = matrix(unlist(df[,-c(1:2)]), nrow=length(russet.map))-1

# simulate data
ploidy=4
np=3
pedigree = diallel_pedigree(parents=np,selfs=0,popsize=200)
pedigreesimR(russet.map, russet.hap[,1:(np*ploidy)],
             ploidy=ploidy,
             sampleHap = FALSE,filename="",
             pedigree,workingfolder = getwd(),
             epsilon = c(0.01, 0.01),
             missingFreq=c(0.1,0.1)
)

# clear up and rename outfiles
# delete all pedsim*.* files. 
sapply(list.files(pattern = "pedsim*"), unlink)
# sapply(list.files(pattern = "GBS*"), unlink)
unlink("*truevalue.csv")
unlink("*geno.csv")
file.rename(list.files(pattern="*geno_snparray.csv")[1], "sim4x_geno.csv")
file.rename(list.files(pattern="*truevalue_ancestral.csv")[1], "sim4x_true.csv")
file.rename(list.files(pattern="*pedigree.csv")[1], "sim4x_ped.csv")
