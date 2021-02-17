##
## Stage II: Variant calling using STITCH
## follows read data processing (readDataProcessing.s)
## followed by Beagle imputation (imputation_beagle.R)
##

##
## split chromosome into pieces to speed up
##
psf<- function(pos, binSize=7e+6, minSnp=1000, quiet=FALSE){
# binSize: bin size
# minSnp: minmum number of SNPs in each bin
   nPt<- floor((max(pos) - min(pos))/binSize)
      #nPt<- round(nPt/nCpu)
      if(nPt < 1) nPt<- 1
      #nPt<- nPt * nCpu
   ps<- seq(min(pos), max(pos), length.out=nPt+1)
      ps<- round(ps)
      ps[1]<- min(pos)-1
      ps[nPt+1]<- max(pos)
   nS<- NULL
   for(n in 1:nPt){
      nS<- c(nS, sum(pos >= ps[n]+1 & pos <= ps[n+1]))
   }
   ms<- min(minSnp, quantile(nS,0.25))
      ms<- ceiling(ms)
   ok<- FALSE; if(nPt < 2) ok<- TRUE
   while(!ok){# merge small ones that are next to each other
      ok<- TRUE
      for(n in 1:(nPt-1)){
         if(nS[n] < minSnp && nS[n+1] < minSnp){
            ps<- ps[-n-1]
            nS[n]<- nS[n] + nS[n+1]
               nS<- nS[-n-1]
            nPt<- nPt - 1

            if(nPt > 1) ok<- FALSE
            break
         }
      }
   }
   ok<- FALSE; if(nPt < 2) ok<- TRUE
   while(!ok){# merge a small one to its neighbor
      ok<- TRUE
      for(n in 1:nPt){
         if(nS[n] < minSnp){
            if(n == 1){
               nS[n+1]<- nS[n] + nS[n+1]
            }else if(n == nPt){
               nS[n-1]<- nS[n] + nS[n-1]
            }else{
               if(nS[n-1] > nS[n+1]){
                  nS[n+1]<- nS[n] + nS[n+1]
               }else{
                  nS[n-1]<- nS[n] + nS[n-1]
               }
            }
            nS<- nS[-n]
            ps<- ps[-n]
            nPt<- nPt - 1

            ok<- FALSE
            break
         }
      }
   }

   if(!quiet){
      cat("Number of jobs: ", nPt, "\n", sep="")
      cat("Bin size (MB):\n"); print(summary(diff(ps)/1e+6))
   }

   list(nPt = nPt, ps = ps)
}

library(parallel) # enable parallel computing in R
library(STITCH)

args <- commandArgs(TRUE) # catch task ID, tid, from script
args
ch <- as.numeric(args[1]) # args[1] to specific chromosome
if(ch == 21){
   ch<- "X"
}else if(ch == 22){
   ch<- "Y"
}else if(is.element(ch, 1:20)){
   ch<- sprintf("%d", ch)
}else{
   stop("Wrong chromsome specified")
}
out_path <- args[2] # args[2] to specific data output path
refPnls <- args[3] # args[3] to specific reference panel path
if(length(args)==4){
   in_path <- args[4] # args[4] to specific data input path
} else {
   in_path <- args[4:length(args)] # args[4:] to specific data input path
}

##
## Arguments in STITCH - modify as appropriate
##
K <- 16
nGen <- 100
buffer<- 1e+6
niterations<- 1

server_environment <- "server"
inputBundleBlockSize <- NA

## create a temporary directory
tempdir <- paste(out_path, "/Tmp", ch, sep="")
datadir <- file.path(tempdir, paste("chr", ch, "Data", sep=""))
resultdir <- file.path(tempdir, paste("chr", ch, "Result", sep=""))
if(!file.exists(datadir))
   system(paste0("mkdir -p ", datadir))
if(!file.exists(resultdir))
   system(paste0("mkdir -p ", resultdir))
system(paste0("rm -r ", resultdir, "/*", sep=""), ignore.stderr = TRUE)

## Choose bam files for STITCH
if(length(in_path)==1){
   bams <- system(paste("ls", paste0(in_path, "/*_mkDup.bam")), intern=TRUE)
} else {
   bams <- unlist(lapply(in_path, function(x) system(paste("ls", paste0(x, "/*_mkDup.bam")), intern=TRUE)))
}

bamlist <- file.path(datadir, paste("bamlist", ch, ".txt", sep=""))
write.table(
   bams,
   file = bamlist,
   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
)
smpl<- sapply(strsplit(bams,"-"), function(x) x[3])
if(sum(table(smpl) != 1) >0) stop("Something wrong")
smplNames <- file.path(datadir, paste("smplNames", ch, ".txt", sep=""))
write.table(
   smpl,
   file = smplNames,
   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
)
genfile <- ""
# ref<- "/projects/ps-palmer/reference_genomes/rat/rn6.fa"

##
## choose founder haplotype data for STITCH
##

refHap<- paste(refPnls, "/phasedVariants_chr", ch, ".hap.gz", sep="")
refLgd<- paste(refPnls, "/phasedVariants_chr", ch, ".legend.gz", sep="")

## physical map
phyMap<- read.table(refLgd,header=TRUE,stringsAsFactors=FALSE)
phyMap<- cbind(sapply(strsplit(phyMap[,"id"],":"),function(x)x[1]),phyMap[,c(2,1,3,4)])
colnames(phyMap)<- c("chr","pos","id","ref","alt")
   #idx<- duplicated(phyMap[,"pos"])
   #phyMap<- phyMap[!idx,]
summary(phyMap[,"pos"])
posfile <- file.path(datadir, paste("pos", ch, ".txt", sep=""))
write.table(
   phyMap[,-3],
   file = posfile,
   sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE
)

ps<- psf(phyMap[,"pos"], binSize=7e+6, minSnp=10000, quiet=FALSE)
   nPt<- ps$nPt
   ps<- ps$ps
nCpu<- min(nPt,5)
nCore<- 1

#msg <- paste(resultdir, "/msg.txt", sep="")
#con<- file(msg, open="w+")
#sink(con, type="message", split=FALSE, append=TRUE)

##
## run STITCH in parallel
##
lst<- list()
nes<- 0
pid<- NULL
cnt<- 0
for(n in 1:nPt){
   cnt<- cnt+1
   pid<- c(pid, mcparallel(
   {###---------------------
      cat("stitch: now process ", n, " out of ", nPt, " blocks. Please wait...\n", sep="")
      if(grepl("rn7", refPnls,fixed=TRUE) && ch == "X"){
         temp_ch <- "X"
      } else if(grepl("rn7", refPnls,fixed=TRUE) && ch == "Y"){
         temp_ch <- "Y"
      }else{
         temp_ch <- paste("chr",ch,sep="")
      }

      tmd<- paste(resultdir, "/tmd", sprintf("%03d", n, sep=""), sep="")
      if(!file.exists(tmd)){
         system(paste0("mkdir -p ", tmd))
      }else{
         system(paste("rm -rf ", tmd, sep=""))
         system(paste0("mkdir -p ", tmd))
      }

      rsd <- file.path(resultdir, paste("rsd", sprintf("%03d",n), sep=""))
      if(!file.exists(rsd)){
         system(paste0("mkdir -p ", rsd))
      }else{
         system(paste("rm -rf ", rsd, sep=""))
         system(paste0("mkdir -p ", rsd))
      }

      st<- STITCH(
                 regionStart = ps[n]+1,
                   regionEnd = ps[n+1],
                      buffer = buffer,
                      method = "diploid",
                   outputdir = rsd,
                         chr = temp_ch,
                     posfile = posfile,
                     genfile = genfile,
                     bamlist = bamlist,
            sampleNames_file = smplNames,
                   #reference = ref, # no need for bams
    reference_haplotype_file = refHap,
      #reference_sample_file = refSmpl,
       reference_legend_file = refLgd,
                           K = K,
                 niterations = niterations,
  shuffleHaplotypeIterations = NA, ## disable heuristics for niterations=1
            refillIterations = NA, ## disable heuristics for niterations=1
                     tempdir = tmd,
                      nCores = nCore,
                        nGen = nGen,
        inputBundleBlockSize = inputBundleBlockSize
      )
      system(paste("/projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools index -f -t ", rsd, "/stitch.", temp_ch,"*.vcf.gz", sep=""))
      cat("\n=>", sprintf("%3d",n), "/", nPt, " is done!\n", sep="")
      system(paste("rm -rf ", tmd, sep=""))
      
      paste("Job #", n, ": okay!", sep="")
   },###---------------------
      mc.set.seed = FALSE,
      silent = FALSE
   )$pid)
   Sys.sleep(5)
   if(n == nPt){
      oTmp<- mccollect(pid,wait=TRUE)
      # process...
      lst<- c(lst, oTmp)
      nes<- sum(nes, length(grep("Error",oTmp)))
      gc()
   }else if(cnt %% nCpu == 0){
      while(TRUE){# wait for some job to finish
         Sys.sleep(30)
         idxTmp<- NULL
         for(mTmp in 1:cnt){
            oTmp<- mccollect(pid[mTmp],wait=FALSE)
            if(!is.null(oTmp[[1]])){
               # process...
               lst<- c(lst, oTmp)
               idxTmp<- c(idxTmp, mTmp)
               nes<- sum(nes, length(grep("Error",oTmp)))
            }
         }
         if(length(idxTmp) > 0){
            cnt<- cnt - length(idxTmp)
            pid<- pid[-idxTmp]
            break
         }
         rm(idxTmp, mTmp, oTmp)
      }
      gc()
   }
}

lst

cat("********Number of failed blocks:", nes, "********\n")

fls<- NULL
for(n in 1:nPt){
   f<- list.files(paste(resultdir, "/rsd", sprintf("%03d",n), sep=""), "stitch.*.vcf.gz$")
   if(length(f) != 1){
      cat(n, " ")
      next
   }
   f<- paste(resultdir, "/rsd", sprintf("%03d",n),"/",f,sep="")
   fls<- c(fls, f)
}
cat("...Move on...\n")
if(length(fls) != nPt) stop("File counts wrong")
str<- "/projects/ps-palmer/software/local/src/bcftools-1.10.2/bcftools concat --no-version -a -d none"
   str<- paste(str, " -O z -o ", out_path, "/chr", ch, "_hs_stitch.vcf.gz", sep="")
   str<- paste(str, " ", paste(fls, collapse=" "), sep="")
system(str, ignore.stdout=FALSE)

str<- paste("rm -rf ", tempdir, sep="")
system(str)

q("no")

##################
# the end #
###########
