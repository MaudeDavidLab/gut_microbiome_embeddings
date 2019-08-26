
library(dada2); packageVersion("dada2")
# Filename parsing
path <- "/nfs3/PHARM/David_Lab/christine/public_data/AG_new/fastqs/runs_use" # CHANGE ME to the directory containing your demultiplexed fastq files

filtpath_100 <- file.path(path, "filtered_100") # Filtered files go into the filtered/ subdirectory
filtpath_150 <- file.path(path, "filtered_150") # Filtered files go into the filtered/ subdirectory
fns <- list.files(path, pattern="fastq.gz") # CHANGE if different file extensions

for(file in fns){
  p = plotQualityProfile(file.path(path,file))
  num_bp = max(p$data$Cycle)
  print(num_bp)
  if(num_bp >= 140){
    # Filtering
    filterAndTrim(file.path(path,fns), file.path(filtpath_150,fns),
                  truncLen=140, maxEE=1, truncQ= 11, rm.phix=TRUE,
                  compress=TRUE, verbose=TRUE, multithread=FALSE)
  }
  if(num_bp < 140){
    # Filtering
    filterAndTrim(file.path(path,fns), file.path(filtpath_100,fns),
                  truncLen=100, maxEE=1, truncQ= 11, rm.phix=TRUE,
                  compress=TRUE, verbose=TRUE, multithread=FALSE)
  }

}
