test_that("shiny", {
  data(dd_list_original, package = "DeeDee")

  library("shiny")

  myapp <- ddedde_app(deedee_obj = DeeDeeLegacy::DeeDeeObject(DeeDeeList = dd_list_original))

  expect_is(myapp, "shiny.appobj")
})
