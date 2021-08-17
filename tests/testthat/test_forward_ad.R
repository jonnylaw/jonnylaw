f <- function(x) {
  5 * x * x + 3 * x + 10
}

test_that(
  "Derivative of 5x^2 + 3x + 10 is 10a + 3",
  hedgehog::forall(
    list(a = hedgehog::gen.c(hedgehog::gen.element(-100:100))),
    function(a) {
      expect_equal(object = f(variable(a))$eps, 
                   expected = 10 * a + 3, 
                   tolerance = 1e-3)
    }
  )
)
