# dada2_processing.R

library(dada2.pipeline)
library(data.table)
library(magrittr)
library(stringr)
library(readxl)

###### Assign/change these variables as needed ######
raw.seq.dirs <- c(
  Run1 = "/dfs/Sharpton_Lab/keaton/Microbiome_methods_FASTQs/Run1",
  Run2 = "/dfs/Sharpton_Lab/keaton/Microbiome_methods_FASTQs/Run2"
)
inDir <- "."
sample.ids.filename <-  "all_attributes_withRunIDs.csv"
file.ids.col <- "SeqID"
run.ids.col <- "RunID"
smpl.ids.col <- "Sample"
ids.tbl.cols <- c(file.ids.col, run.ids.col, smpl.ids.col)
split <- "--"
###### Begin script ######
orig.dir <- getwd()
fastq.dir <- file.path(inDir, "FastQs")
if (!dir.exists(fastq.dir)) { dir.create(fastq.dir) }
if (!dir.exists(inDir)) { dir.create(inDir) }

setwd(orig.dir)
sample.ids.dt <- read.csv(file.path(inDir, sample.ids.filename))
if (!("data.table" %in% sample.ids.dt)) {
  sample.ids.dt <- as.data.table(sample.ids.dt) %>%
    setkeyv(smpl.ids.col)
}
###### Symlink FASTQs ######
setwd(fastq.dir)
symlink.fastqs(
  seq.dirs = raw.seq.dirs,
  ids.tbl = sample.ids.dt[, ..ids.tbl.cols],
  smpl.id.col = smpl.ids.col,
  file.id.col = file.ids.col,
  run.id.col = run.ids.col,
  split.pattern = split,
  quiet = FALSE
)
