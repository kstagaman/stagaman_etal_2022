# gen_clr_matrix.R

gen.clr.matrix <- function(
  asv.mat,
  min_reads,
  min_prop = 0.001,
  min_occur = 0,
  smpls_by_row = TRUE,
  method = "CZM",
  lab = 0
  ) {
  require(CoDaSeq)
  require(zCompositions)
  asv.mat.f <- codaSeq.filter(
    asv.mat,
    min.reads = min_reads,
    min.prop = min_prop,
    min.occurrence = min_occur,
    samples.by.row = smpls_by_row
  )
  # replace 0 values with an estimate
  asv.mat.fn0 <- cmultRepl(t(asv.mat.f), method = method, label = lab)
  return(codaSeq.clr(asv.mat.fn0))
}
