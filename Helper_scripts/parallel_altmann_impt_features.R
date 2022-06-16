# parallel_altmann_impt_features.R

parallel.altmann.impt.features <- function(randForest, dependent.var, num.cores, num.permutations = 100) {
  x <- randForest$finalModel
  data <- as.data.frame(randForest$trainingData)
  dat_x <- data[, -c(which(names(data) == dependent.var))]
  dat_x <- dat_x[, names(x$variable.importance)] %>% as.matrix()
  set.seed(42)
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  vimp <- foreach(
    i = 1:num.permutations,
    .final = rlist::list.cbind,
    .verbose = TRUE
  ) %dopar% {
    dat_y <- data[sample(nrow(data)), dependent.var]
    if (class(dat_y) == "character") {dat_y <- factor(dat_y)}
    ranger::ranger(
      x = dat_x,
      y = dat_y,
      num.trees = x$num.trees,
      mtry = x$mtry,
      min.node.size = x$min.node.size,
      importance = x$importance.mode,
      replace = x$replace,
      num.threads = 1
    )$variable.importance
  }
  stopCluster(cl)
  pval <- sapply(1:nrow(vimp), function(i) {
    (sum(vimp[i, ] >= x$variable.importance[i]) + 1)/(ncol(vimp) + 1)
  })

  impt.mat <- cbind(x$variable.importance, pval)
  colnames(impt.mat) <- c("Importance", "Pvalue")

  impt.dt <- as.data.table(impt.mat, keep.rownames = "ASV")
  impt.dt[, Covar := dependent.var]
  return(impt.dt)
}
