dde <- DeeDeeExperiment(
  se_macrophage_noassays,
  de_results = del
)

# needed for bad handling of namespace by ggvenn
library("dplyr")

test_that("plotting & plotting", {
  p_scatter <- deedee_scatter(dde)
  expect_is(p_scatter, "gg")


  p_heatmap <- deedee_heatmap(dde)
  expect_is(p_heatmap, "Heatmap")


  p_venn <- deedee_venn(dde)
  expect_is(p_venn, "gg")


  p_upset <- deedee_upset(dde)
  expect_is(p_upset, "gg")


  p_cat <- deedee_cat(dde)
  expect_is(p_cat, "gg")


  p_qq <- deedee_qq(dde)
  expect_is(p_qq, "gg")
  p_qq_m <- deedee_qqmult(dde)
  expect_is(p_qq_m, "gg")


  # deedee_summary(dde)
})


