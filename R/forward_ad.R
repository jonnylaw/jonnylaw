#' Representation of a Dual Number
#'
#' @param real
#' @param eps
#'
#' @return
#' @export
#'
#' @examples
dual <- function(real, eps) {
  if (!is.numeric(real)) {
    stop("real must be numeric")
  }
  structure(list(real = real, eps = eps), class = "dual")
}

#' Build a variable
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
variable <- function(x) {
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  dual(x, 1)
}

#' Build a constant dual number
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
const <- function(x) {
  if (!is.numeric(x)) {
    stop("x must be numeric")
  }
  dual(x, 0)
}

#' Addition in the dual numbers
#'
#' @param x 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
plus <- function(x, y) {
  dual(x$real + y$real, x$eps + y$eps)
}

#' Minus in the dual numbers
#'
#' @param x 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
minus <- function(x, y) {
  dual(x$real - y$real, x$eps - y$eps)
}


#' Multiplication in the dual numbers
#'
#' @param x 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
times <- function(x, y) {
  dual(x$real * y$real, x$eps * y$real + y$eps * x$real)
}


#' Divide in the dual numbers
#'
#' @param x 
#' @param y 
#'
#' @return
#' @export
#'
#' @examples
divide <- function(x, y) {
  dual(
    x$real / y$real,
    (x$eps * y$real - x$real * y$eps) / (y$real * y$real)
  )
}

#' Lift a function to operate on bare numbers
#'
#'
#' @param f a function operating on dual numbers
#'
#' @return
#'
#' @examples
lift_function <- function(f) {
  function(x, y) {
    if (is.double(x)) {
      f(const(x), y)
    } else if (is.double(y)) {
      f(x, const(y))
    } else {
      f(x, y)
    }
  }
}

#' Operators for univariate Dual Numbers
#'
#' @param e1
#' @param e2
#'
#' @return
#' @export
#'
#' @examples
Ops.dual <- function(x, y) {
  switch(
    .Generic,
    `+` = lift_function(plus)(x, y),
    `-` = lift_function(minus)(x, y),
    `*` = lift_function(times)(x, y),
    `/` = lift_function(divide)(x, y)
  )
}
