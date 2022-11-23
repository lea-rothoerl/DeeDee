test_that("creating", {
  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = del
  )
  dde

  expect_is(dde, "DeeDeeExperiment")

  dde_only_de <- DeeDeeExperiment(
    de_results = del
  )
  expect_is(dde, "DeeDeeExperiment")

})


test_that("adding and removing", {
  dde <- DeeDeeExperiment(
    se_macrophage_noassays,
    de_results = del
  )

  new_del <- list(
    ifng2 = del$ifng_vs_naive,
    ifngsalmo2 = del$ifngsalmo_vs_naive
  )
  # add a new (set of) DE result(s)
  dde_new <- add_dea(dde, new_del)
  expect_is(dde_new, "DeeDeeExperiment")
  expect_equal(length(dea(dde)), 4)
  expect_equal(length(dea(dde_new)), 6)

  dde_removed <- remove_dea(dde, "ifngsalmo_vs_naive")
  expect_is(dde_removed, "DeeDeeExperiment")
  expect_equal(length(dea(dde_removed)), 3)
})
