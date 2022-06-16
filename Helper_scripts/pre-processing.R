# pre-processing.R

### LOCAL ###

## BELOW I NEED TO MAKE SURE THAT PLAIN BLK AND WATER SAMPLES (NO KIT IDS) ARE GRABBING THE CORRECT SEQUENCING FILES
## also need to keep sequencing barcodes for seq file IDs
library(data.table)
library(magrittr)
library(stringr)
seq.smpl.ids <- read.csv("Input/seq_sample_ids.csv") %>% as.data.table() %>% unique()
seq.smpl.ids[, RunID := sapply(str_split(SeqID, "-"), tail, 1)]
seq.smpl.ids[, SeqID := str_remove(SeqID, "-R[12]")]
seq.smpl.ids[, SeqID := str_replace_all(SeqID, "-", "_")]
seq.smpl.ids[, SeqID := str_replace_all(SeqID, " ", "-")]
seq.smpl.ids[, Sample := SeqID]
seq.smpl.ids[, Sample := ifelse(str_detect(Sample, "^[0-9]"), paste0("P", Sample), Sample)]
seq.smpl.ids[, Sample := str_replace(Sample, "([A-H])([0-9][\\._])", "\\10\\2")]
seq.smpl.ids[, OriginID := sapply(str_split(Sample, "_"), extract, 1)]
dissect.res <- read.csv("Input/dissection_results.csv") %>% as.data.table()
dup.wells <- str_subset(seq.smpl.ids$OriginID, "\\.") %>%
  str_split("\\.") %>%
  sapply(`[`, 1) %>%
  unique()
for (well in dup.wells) {
  row.idx <- which(dissect.res$Well == well)
  row1 <- dissect.res[row.idx]
  row2 <- copy(row1)
  row1$Well <- paste0(row1$Well, ".1")
  row2$Well <- paste0(row2$Well, ".2")
  dissect.res <- dissect.res[-row.idx]
  sid.row.idx <- which(seq.smpl.ids$OriginID == well)
  dissect.res <- rbind(dissect.res, row1, row2)
}
sample.data <- dissect.res[Alive.Dead == "alive"][seq.smpl.ids, on = "Well==OriginID", nomatch = 0]
names(sample.data)[which(names(sample.data) == "Well")] <- "OriginID"
setkey(sample.data, Sample)
sample.data[, Alive.Dead := NULL]
sample.data[, Extraction.Batch := NULL]
sample.data[
  ,
  Extraction.Kit := ifelse(
    is.na(Extraction.Kit) & str_detect(OriginID, "-"),
    ifelse(str_detect(OriginID, "BT"), "Blood & Tissue", ifelse(str_detect(OriginID, "MN"), "NucleoSpin", "PowerSoil")),
    ifelse(Extraction.Kit == "Kit1", "PowerSoil", ifelse(Extraction.Kit == "Kit2", "NucleoSpin", "Blood & Tissue"))
  )
]
sample.data[
  ,
  Single.Triplicate := ifelse(str_detect(Sample, "_S_"), "PCRx1", ifelse(str_detect(Sample, "_T_"), "PCRx3", NA))
  ]
sample.data[, Single.Type := ifelse(str_detect(Sample, "_1$"), "subset", "separate")]
sample.data[, SeqID := str_remove(SeqID, "-")]
sample.data[, SeqID := str_replace_all(SeqID, "_", "-")]
# sample.data <- unique(sample.data)
saveRDS(sample.data, file = "Input/dt_sample_data.rds")

#################################################################################

### DARWIN ###

cleanup <- FALSE
test.run <- FALSE
use.generated.output.path <- TRUE
alt.output.path <- ""
test.run.nSamples <- 100
cores <- 50
raw.seq.dir1 <- "/nfs2/hts/miseq/210205_M01498_0763_000000000-JG3KB/Data/Intensities/BaseCalls"
raw.seq.dir2 <- "/nfs2/hts/miseq/210212_M01498_0765_000000000-JG37D/Data/Intensities/BaseCalls"
raw.seq.dirs <- c(raw.seq.dir1, raw.seq.dir2)
split.pattern <- "--"

source("Helper_scripts/dada2_processing_functions.R")
setDTthreads(threads = cores)
sample.data <- readRDS("Input/dt_sample_data.rds")
orig.dir <- getwd()
fastq.dir <- "Input/FastQs"
if (!dir.exists(fastq.dir)) {dir.create(fastq.dir)}
if (use.generated.output.path) {
  output <- output.path
} else {
  output <- alt.output.path
}

setwd(fastq.dir)
seqIDs <- sample.data$SeqID
if (test.run) {
  set.seed(42)
  seqIDs <- sample(
    sample.data[RunID == "R1" & !str_detect(OriginID, "BLK|Water|^G05")]$SeqID, ###
    test.run.nSamples
    )
}
# seqID <- "2H4-T-1"
for (seqID in seqIDs) {
  seq.dir <- ifelse(sample.data[SeqID == seqID]$RunID == "R1", raw.seq.dirs[1], raw.seq.dirs[2])
  barcode <- sample.data[SeqID == seqID]$Barcode
  cat(seqID, sep = "\n")
  seq.files <- list.files(
    path = seq.dir,
    pattern = barcode,
    full.names = T
  )
  if (length(seq.files) == 0) {stop("no files matched barcode")}
  if (test.run) {
    while (length(seq.files) > 2) {
      newID <- sample(sample.data[RunID == "R1" & !(SeqID %in% seqIDs)]$SeqID, 1)
      seq.files <- list.files( ### alter matching to barcode seqs
        path = raw.seq.dir1,
        pattern = newID,
        full.names = T
      )
    }
  } else {
    if (length(seq.files) > 2) {
      stop(paste0(seqID, ": Too many matching files!"))
    }
  }
  read1.file <- seq.files[str_detect(seq.files, "_R1_")]
  read2.file <- seq.files[str_detect(seq.files, "_R2_")]
  lnName.r1 <- paste(
    sample.data[SeqID == seqID]$Sample,
    "R1.fastq.gz",
    sep = split.pattern
  )
  lnName.r2 <- paste(
    sample.data[SeqID == seqID]$Sample,
    "R2.fastq.gz",
    sep = split.pattern
  )
  cmd1 <- paste("ln -s", read1.file, lnName.r1)
  cmd2 <- paste("ln -s", read2.file, lnName.r2)
  system(cmd1)
  system(cmd2)
}
setwd(orig.dir)

dada2.upto.qualPlots(
  fastq.path = fastq.dir,
  file.split.pattern = "--",
  maxCores = cores,
  random.seed = 42
)

system(paste0("cp -v dada2*", Sys.Date(), "_output/*.pdf ~/temp"))
system("cp -v /home/micro/stagamak/databases/Tree-building_files/* ./")

dada2.finish(
  fastq.path = fastq.dir,
  truncFwd.len = 250,
  truncRev.len = 200,
  taxa.db.name = "silva_nr_v138",
  metadata.file = "Input/dt_sample_data.rds",
  guide.seqs.file = "100_silva_123_pd.fasta",
  alignment.template.file = "silva.seed_v138.align",
  maxCores = cores,
  build.tree = TRUE
)

system(paste("mv mothur*", output))
if (file.exists("tmp.txt")) {file.remove("tmp.txt")}
file.copy(
  from = file.path(output, "phyloseq.rds"),
  to = "Input/phyloseq.rds",
  overwrite = TRUE
)

if (cleanup) {
  system(paste0("tar zvcf ", output, ".tgz ", output))
  system(paste("rm -r", output))
  system(paste("rm -r", fastq.dir))
  system("rm *silva*")
}
