# dada2_processing_functions.R

#### TO DO: REDIRECT ALL OUTPUT TO THE DADA2 DIRECTORY
library(dada2)
library(ggplot2)
library(cowplot)
library(phyloseq)
library(phyloseqCompanion)
library(seqinr)
library(stringr)
library(phangorn)
library(ape)
theme_set(theme_cowplot())

getN <- function(x) sum(getUniques(x))
my.cat <- function(x) {cat(paste0("\n### ", x, "\n"), sep = "\n")}

dada.pkg.ver <- paste0("dada2_", packageVersion("dada2"))
process.date <- Sys.Date()
process.id <- paste(dada.pkg.ver, process.date, sep = "_")
output.path <- paste0(process.id, "_output")
if (!dir.exists(output.path)) {
  dir.create(output.path)
}
proj.name <- tail(strsplit(getwd(), "/")[[1]], 1)

run.env <- new.env()

symlink.fastqs <- function(seq.dirs, link.dir, ids.dt, split.pattern = "_") {
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
}

# this function runs up through making quality plots, allowing the user to manually choose the truncation lengths based on those results

dada2.upto.qualPlots <- function(
  fastq.path,
  file.split.pattern = "-",
  maxCores = 4,
  random.seed = 42,
  user.output.path = NULL,
  force = FALSE
  ) {
  print(head(list.files(fastq.path), 20))
  proceed <- readline(
    prompt = "Above are the first 20 (or fewer) files in the provided path, do you want to proceed? [y/n]: "
  )
  while (!(proceed %in% c("y", "n"))) {
    proceed <- readline(prompt="Please answer y or n: ")
  }
  if (proceed == "n") {
    my.cat("Terminated")
  } else {
    if (is.null(user.output.path)) {
      output <- output.path
    } else {
      output <- user.output.path
    }
    run.env$fnFs <- list.files(
      fastq.path,
      pattern = paste0(file.split.pattern, "R1"),
      full.names = TRUE
    ) %>% sort()
    run.env$fnRs <- list.files(
      fastq.path,
      pattern = paste0(file.split.pattern, "R2"),
      full.names = TRUE
    ) %>% sort()
    run.env$sample.names <- strsplit(
      basename(run.env$fnFs),
      file.split.pattern
      ) %>% sapply(`[`, 1)
    if (length(run.env$sample.names) > 20) {
      set.seed(random.seed)
      plot.samples <- sort(sample(run.env$sample.names, 20, replace = F))
      plot.fnFs <- run.env$fnFs[
        str_detect(run.env$fnFs, paste(plot.samples, collapse = "|"))
      ]
      plot.fnRs <- run.env$fnRs[
        str_detect(run.env$fnRs, paste(plot.samples, collapse = "|"))
      ]
      plot.fnFs.file <- "fwdReads_subset20_qualPlot.pdf"
      plot.fnRs.file <- "revReads_subset20_qualPlot.pdf"
    } else {
      plot.fnFs <- run.env$fnFs
      plot.fnRs <- run.env$fnRs
      plot.fnFs.file <- "fwdReads_qualPlot.pdf"
      plot.fnRs.file <- "revReads_qualPlot.pdf"
    }

    plot.files <- c(file.path(output, plot.fnFs.file), file.path(output, plot.fnRs.file))
    if (all(file.exists(plot.files)) & !force) {
      my.cat("Quality plots already exist, skipping...")
    } else {
      my.cat("Making quality plots...")
      fnFs.qualPlot0 <- plotQualityProfile(plot.fnFs)
      fnFs.qualPlot <- fnFs.qualPlot0 +
        scale_x_continuous(breaks = seq(0, 300, 25)) +
        background_grid(major = "x", minor = "x") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      ggsave(
        fnFs.qualPlot,
        file = plot.files[1]
      )

      fnRs.qualPlot0 <- plotQualityProfile(plot.fnRs)
      fnRs.qualPlot <- fnRs.qualPlot0 +
        scale_x_continuous(breaks = seq(0, 300, 25)) +
        background_grid(major = "x", minor = "x") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
      ggsave(
        fnRs.qualPlot,
        file = plot.files[2]
      )
    }
    my.cat("\tDONE")
  }
}

dada2.finish <- function(
  fastq.path,
  truncFwd.len,
  truncRev.len,
  taxa.db.name,
  metadata.file,
  guide.seqs.file,
  alignment.template.file,
  maxCores = 4,
  build.tree = TRUE,
  user.output.path = NULL,
  force = FALSE
) {
  if (length(names(run.env)) == 0) {
    stop("Fuction 'dada2.upto.qualPlots()' must be run first.")
  }
  if (build.tree) {
    file.check <- list.files(pattern = "silva")
    if (!(
      guide.seqs.file %in% file.check &
      alignment.template.file %in% file.check
    )) {
      stop(
        paste(
          "Files", guide.seqs.file, "and", alignment.template.file,
          "must be in your current directory to build a tree."
        )
      )
    }
  }
  if (is.null(user.output.path)) {
    output <- output.path
  } else {
    output <- user.output.path
  }
  filtFs <- file.path(
    fastq.path,
    "Filtered",
    paste0(run.env$sample.names, "_F_filt.fastq.gz")
  )
  filtRs <- file.path(
    fastq.path,
    "Filtered",
    paste0(run.env$sample.names, "_R_filt.fastq.gz")
  )
  names(filtFs) <- run.env$sample.names
  names(filtRs) <- run.env$sample.names

  out.file <- file.path(output, "filter_and_trim_numbers.rds")
  if (file.exists(out.file) & !force) {
    my.cat("Filtering and trimming completed previously, skipping...")
    out <- readRDS(out.file)
    my.cat("\tDONE:")
    print(head(out))
  } else {
    my.cat("Filtering and trimming...")
    out <- filterAndTrim(
      run.env$fnFs, filtFs,
      run.env$fnRs, filtRs,
      truncLen = c(truncFwd.len, truncRev.len),
      maxN = 0,
      maxEE = c(2, 2),
      truncQ = 2,
      rm.phix = TRUE,
      compress = TRUE,
      multithread = maxCores
    )
    my.cat("\tDONE:")
    print(head(out))
    saveRDS(out, file = out.file)
  }


  filtFs <- filtFs[file.exists(filtFs)]
  filtRs <- filtRs[file.exists(filtRs)]

  err.files <- c(file.path(output, "errF.rds"), file.path(output, "errR.rds"))
  err.plot.files <- c(file.path(output, "errF_plot.pdf"), file.path(output, "errR_plot.pdf"))

  if (all(file.exists(err.files)) & !force) {
    my.cat("Errors learned previously, skipping...")
    errF <- readRDS(err.files[1])
    errR <- readRDS(err.files[2])
  } else {
    my.cat("Learning errors and making error plots...")
    errF <- learnErrors(filtFs, multithread = maxCores)
    saveRDS(errF, file = err.files[1])
    errF.plot <- plotErrors(errF, nominalQ = TRUE)
    ggsave(errF.plot, file = err.plot.files[1])

    errR <- learnErrors(filtFs, multithread = maxCores)
    saveRDS(errR, file = err.files[2])
    errR.plot <- plotErrors(errR, nominalQ = TRUE)
    ggsave(errR.plot, file = err.plot.files[2])
  }
  my.cat("\tDONE")

  dada.files <- c(file.path(output, "dadaFs.rds"), file.path(output, "dadaRs.rds"))
  if (all(file.exists(dada.files)) & !force) {
    my.cat("dada-ing done previously, skipping...")
    dadaFs <- readRDS(dada.files[1])
    dadaRs <- readRDS(dada.files[2])
  } else {
    my.cat("dada-ing...")
    dadaFs <- dada(filtFs, err = errF, multithread = maxCores)
    saveRDS(dadaFs, file = dada.files[1])
    dadaRs <- dada(filtRs, err = errR, multithread = maxCores)
    saveRDS(dadaRs, file = dada.files[2])
  }
  my.cat("\tDONE")

  seqtab.file <- file.path(output, "seqtab.rds")
  mergers.file <- file.path(output, "mergers.rds")

  if (file.exists(seqtab.file) & !force) {
    my.cat("Merging pairs completed previously, skipping...")
    seqtab <- readRDS(seqtab.file)
  } else {
    my.cat("Merging pairs...")
    mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
    merge.results <- sapply(mergers, getN)
    my.cat("\tDONE")

    if (sum(merge.results == 0) >= 3) {
      my.cat("Multiple sequence merging fails, proceeding with FORWARD READS ONLY")
      seqtab <- makeSequenceTable(dadaFs)
      saveRDS(seqtab, file = seqtab.file)
    } else {
      saveRDS(mergers, file = mergers.file)
      seqtab <- makeSequenceTable(mergers)
      saveRDS(seqtab, file = seqtab.file)
    }
  }

  seqtab.nochim.file <- file.path(output, "seqtab_nochim.rds")

  if (file.exists(seqtab.nochim.file) & !force) {
    my.cat("Bimeras previously removed, skipping...")
    seqtab.nochim <- readRDS(seqtab.nochim.file)
    my.cat(paste("\t\tPercent seqs kept:", sum(seqtab.nochim)/sum(seqtab)))
  } else {
    my.cat("Removing bimeras...")
    seqtab.nochim <- removeBimeraDenovo(
      seqtab,
      method = "consensus",
      multithread = maxCores,
      verbose = TRUE
    )
    my.cat("\tDONE")
    my.cat(paste("\t\tPercent seqs kept:", sum(seqtab.nochim)/sum(seqtab)))
    saveRDS(
      seqtab.nochim,
      file = seqtab.nochim.file
    )
  }
  if (sum(merge.results == 0) >= 3) {
    track <- cbind(
      out,
      sapply(dadaFs, getN),
      sapply(dadaRs, getN),
      rowSums(seqtab.nochim)
    )

    colnames(track) <- c(
      "input",
      "filtered",
      "denoisedF",
      "denoisedR",
      "nonchimF"
    )
  } else {
    track <- cbind(
      out,
      sapply(dadaFs, getN),
      sapply(dadaRs, getN),
      merge.results,
      rowSums(seqtab.nochim)
    )

    colnames(track) <- c(
      "input",
      "filtered",
      "denoisedF",
      "denoisedR",
      "merged",
      "nonchim"
    )
  }
  rownames(track) <- run.env$sample.names
  my.cat("See 'track.rds' in output directory for how many seqs made it through. First 6 samples look as follows:")
  print(head(track))
  saveRDS(track, file = file.path(output, "track.rds"))

  taxa.file <-  file.path(output, paste0("taxa_", taxa.db.name, ".rds"))

  if (file.exists(taxa.file) & !force) {
    my.cat("Taxonomy previously assigned, skipping...")
    taxa <- readRDS(taxa.file)
  } else {
    my.cat("Assigning taxonomy...")
    taxa <- assignTaxonomy(
      seqtab.nochim,
      paste0("/home/micro/stagamak/databases/", taxa.db.name, "_train_set.fa.gz"),
      multithread = maxCores
    )
    saveRDS(taxa, file = taxa.file)
  }

  my.cat("\tDONE")
  if (str_detect(metadata.file, "\\.csv$")) {
    smpl.df <- read.csv(metadata.file, row.names = 1)
  } else if (str_detect(metadata.file, "\\.tsv$|\\.txt$")) {
    smpl.df <- read.table(
      metadata.file,
      sep = "\t",
      header = TRUE,
      row.names = 1
      )
  } else if (str_detect(metadata.file, "\\.rds$")) {
    smpl.tbl <- readRDS(metadata.file)
    smpl.col <- smpl.tbl[
      ,
      sapply(
        .SD,
        function(col) {
          sum(col %in% run.env$sample.names) == length(run.env$sample.names)
        }
      )
    ] %>%
      extract(., .) %>%
      names()
    if ("data.table" %in% class(smpl.tbl)) {
      smpl.df <- as.data.frame(smpl.tbl)
      row.names(smpl.df) <- smpl.tbl[[smpl.col]]
    } else {
      smpl.df <- smpl.tbl
    }
  } else {
    stop("Meta-data filetype not recognizable from extension")
  }
  smpl.df <- smpl.df[run.env$sample.names, ]
  ps0 <- phyloseq(
    otu_table(seqtab.nochim, taxa_are_rows = FALSE),
    sample_data(smpl.df),
    tax_table(taxa)
  )

  ps1 <- numbered.ASVs(
    ps = ps0,
    # prefix = paste0(proj.name, "_ASV"),
    save.dir = output,
    save.file = "asv_sequences"
  )
  asv.seqs <- readRDS(file.path(output, "asv_sequences.rds"))

  seqinr::write.fasta(
    sequences = as.list(asv.seqs),
    names = taxa_names(ps1),
    as.string = TRUE,
    file.out = file.path(output, "asv_sequences.fasta")
  )
  if (build.tree) {
    my.cat("Proceeding with phylogenetic tree:")
    asv.seqs.file <- file.path(output, "asv_sequences.fasta")
    asv.withguides.file  <- file.path(output, "asv_and_guide_seqs.fasta")
    asv.tree.rooted.file <- file.path(output, "asv_NASTaligned_seqs.nwk")

    cmd <- paste(
      "cat",
      asv.seqs.file,
      guide.seqs.file,
      ">",
      asv.withguides.file
    )
    system(cmd)

    my.cat("Aligning sequences...")
    cmd <- paste0(
      "mothur \"#align.seqs( fasta=",
      asv.withguides.file,
      ", reference=",
      alignment.template.file,
      ", flip=t",
      ", keepdots=t",
      ", processors=", maxCores,
      ", outputdir=",
      output,
      "/ )\""
    )
    system(cmd)
    my.cat("\tDONE")

    mothur.output.file <- file.path(output, "asv_and_guide_seqs.align")
    fasttree.log.file <- file.path(output, "fasttree.log")
    fasttree.output.file <- file.path(output, "asv_and_guide_seqs.nwk")

    my.cat("Building phylogenetic tree...")
    cmd <- paste0(
      "export OMP_NUM_THREADS=",
      maxCores,
      "; $HOME/src/FastTreePar-2.1.9 -nt -nosupport -quote -gtr -gamma -log ",
      fasttree.log.file,
      " ",
      mothur.output.file,
      " > ",
      fasttree.output.file

    )
    system(cmd)
    my.cat("\tDONE")
    asvs.and.guides.tree <- read_tree(fasttree.output.file)
    asvs.and.guides.tree.rooted <- phangorn::midpoint(asvs.and.guides.tree)

    guides <- scan(guide.seqs.file, what = "character" )
    guide.ids <- guides[stringr::str_detect(guides, ">" )]
    guide.ids <- stringr::str_remove(guide.ids, ">")

    asvs.tree.rooted <- ape::drop.tip(asvs.and.guides.tree.rooted, guide.ids)
    write.tree(asvs.tree.rooted, file = asv.tree.rooted.file)

    phy_tree(ps1) <- phy_tree(asvs.tree.rooted)
  }
  saveRDS(ps1, file = file.path(output, "phyloseq.rds"))
  my.cat("\tDONE and DONE")
}
