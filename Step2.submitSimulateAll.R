###############################################
## pipeline to submit simulation to cluster for all simulations
## for SimClone-1000 simulations
## author: maxime.tarabichi@crick.ac.uk, 2017-2018
## for PCAWG-11
###############################################


###############################################
## indices of simulations
IDS <-  1:700
###############################################



###############################################
######### create slurm submit file for one sample
generateSubmitFile<-function(sampIndex,file)
{
    cat("#!/bin/bash",file=file)
    cat("\n",file=file,append=T)
    cat("## a cluster job for simclone",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("#SBATCH -n 1",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("#SBATCH -t 5-24:00:00 ",file=file,append=T)## this is too much for most, but some samples take a long time to simulate, we removed them if taking more than 2 days
    cat("\n",file=file,append=T)
    cat(paste("#SBATCH --job-name=sim",sampIndex,sep=""),file=file,append=T)
    cat("\n",file=file,append=T)
    cat("#SBATCH --mem=15G",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("\n",file=file,append=T)
    cat("## execute the simclone pipeline",file=file,append=T)
    cat("\n",file=file,append=T)
    cat(paste("module load R"),
        file=file,append=T)
    cat("\n",file=file,append=T)
    cat(paste("module load R-bundle-Bioconductor/3.2-foss-2016a-R-3.2.3"),
        file=file,append=T)
    cat("\n",file=file,append=T)
    cat(paste("Rscript simulateSample.R",
	      sampIndex, #1
              sep=" "),
        file=file,append=T)
    cat("\n",file=file,append=T)
    system(paste("chmod u+rwx ",file,sep=""))
}

callBash<-function(file)
{
    print(file)
    print(system(paste("sbatch ",file," -e ~/err.txt &",sep="")))
}


submitJob <- function(id)
{
    submitFile <- paste0("~/submitFiles/",id,".final.bash")
    generateSubmitFile(sampIndex=id,
                       file=submitFile)
    callBash(submitFile)
}
###############################################


###############################################
## submit ALL
kNull <- lapply(IDS,submitJob)
###############################################
q(save="no")
###############################################



