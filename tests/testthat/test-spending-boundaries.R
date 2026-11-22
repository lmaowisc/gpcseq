test_that("four spending families reproduce first-crossing probabilities at correlated looks", {
  # Independent normal integration: Z_j = l_j X + sqrt(1-l_j^2) E_j,
  # where X and E_j are independent standard normals. Integrating the product
  # of conditional continuation probabilities over X yields these boundaries.
  loading <- c(.725, .85)
  correlation <- outer(loading, loading)
  diag(correlation) <- 1
  families <- c("obrien-fleming", "pocock", "power", "hsd")
  final <- c(1.98676131835109, 2.23985511843361, 2.03399558724926, 1.98905783818126)
  for (j in seq_along(families)) {
    plan <- if (families[j] == "power")
      gpc_design(c(1, 2), spending = "power", spending_param = 2)
    else gpc_design(c(1, 2), spending = families[j])
    planned <- plan$design$spending
    c1 <- .gpc_boundary(diag(1), numeric(), "two-sided", planned[1])$critical
    c2 <- .gpc_boundary(correlation, c1, rep("two-sided", 2), planned[2])$critical
    expect_equal(c2, final[j], tolerance = 2e-6)
  }
})
