# random_forest_functions.R

gen.bap.random.forest <- function(physeq, rnd.seed) {
  feature.dt <- as.data.table(otu.matrix(physeq), keep.rownames = "Sample") %>%
    setkeyv("Sample")
  bap.dt <- sample.data.table(physeq)[, .(Sample, BaP.Treat)] %>% setkeyv("Sample")
  rf.data <- bap.dt[feature.dt]
  rf.data[, Sample := NULL]

  bap.recipe <- recipe(rf.data) %>%
    update_role(BaP.Treat, new_role = "outcome") %>%
    update_role(-BaP.Treat) %>%
    step_nzv(all_predictors())

  set.seed(rnd.seed)
  cl <- makeCluster(num.cores, type = "FORK", outfile = "")
  registerDoParallel(cl, num.cores)
  mod <- train(
    bap.recipe,
    data = as.data.frame(rf.data),
    method = "ranger",
    importance = "permutation",
    metric = "ROC",
    tuneGrid = expand.grid(
      mtry = c(10, 25, 50, 100),
      splitrule = c("gini", "extratrees"),
      min.node.size = 1 # 1 for classification, 5 for regression according to `ranger` manual
    ),
    trControl = trainControl(
      method = "cv",
      number = 5,
      summaryFunction = twoClassSummary, # to calculate ROC
      classProbs = TRUE, # IMPORTANT for twoClassSummary
      savePredictions = "all",
      returnResamp = "all",
      verboseIter = TRUE
    )
  )
  stopCluster(cl)
  return(mod)
}

