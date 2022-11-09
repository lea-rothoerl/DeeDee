dde <- DeeDeeExperiment(
  se_macrophage_noassays,
  de_results = del
)

test_that("plotting & plotting", {
  p_scatter <- ddedde_scatter(dde)
  expect_is(p_scatter, "gg")


  p_heatmap <- ddedde_heatmap(dde)
  expect_is(p_heatmap, "Heatmap")


  p_venn <- ddedde_venn(dde)
  expect_is(p_venn, "gg")


  p_upset <- ddedde_upset(dde)
  expect_is(p_upset, "gg")


  p_cat <- ddedde_cat(dde)
  expect_is(p_cat, "gg")


  p_qq <- ddedde_qq(dde)
  expect_is(p_qq, "gg")
  p_qq_m <- ddedde_qqmult(dde)
  expect_is(p_qq_m, "gg")


  # ddedde_summary(dde)
})


