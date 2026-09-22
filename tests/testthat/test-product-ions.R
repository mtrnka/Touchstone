test_that("product ions map to unique backbone positions", {
  result <- getProductIonMatches(
    "b2;b3;y2;y3;b4-H2O;unrelated",
    pep.len = 8
  )

  expect_equal(result$ion_indicies, c(2, 3, 5, 6))
  expect_equal(result$n_bonds_possible, 7)
  expect_equal(result$n_bonds_observed, 4)
  expect_equal(result$longest_ladder, 2)
  expect_equal(result$longest_gapped_ladder_observed, 4)
  expect_equal(result$Perc.Gapped.Ladder, 4 / 7)
})

test_that("missing product-ion annotations produce zero evidence", {
  result <- getProductIonMatches(NA_character_, pep.len = 8)

  expect_length(result$ion_indicies, 0)
  expect_equal(result$n_bonds_possible, 7)
  expect_equal(result$n_bonds_observed, 0)
  expect_equal(result$longest_gapped_ladder_observed, 0)
  expect_equal(result$Perc.Gapped.Ladder, 0)
})

test_that("product-ion features support transparent filtering", {
  input <- data.frame(
    MSMS.Ions.1 = c("b1;b2;b3", "b1"),
    MSMS.Ions.2 = c("y1;y2;y3", "y1;y2;y3"),
    Len.Pep.1 = c(5, 5),
    Len.Pep.2 = c(5, 5),
    row_id = c("pass", "fail")
  )

  annotated <- calculateProductIons(input)
  filtered <- productIonFilter(
    annotated,
    minProducts.1 = 3,
    minProducts.2 = 3
  )

  expect_equal(annotated$numProdIons.1, c(3, 1))
  expect_equal(annotated$numProdIons.2, c(3, 3))
  expect_identical(filtered$row_id, "pass")
})

