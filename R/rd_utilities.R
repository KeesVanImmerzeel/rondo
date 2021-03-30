
#' Internal function to determine the number of areas for which kD- and c-values are specified.
#'
#' @param kD Transmissivity (m2/d) (numeric vector).
#' @return The number of areas for which kD- and c-values are specified (-) (integer).
#' @details The number of circular areas is equal to the length of the vector kD minus 1.
.create_m <- function(kD) {
  return(length(kD))
}

#' Internal function to calculate the leakage factor.
#'
#' @inheritParams .create_m
#' @param c Hydraulic resistance of top layer (d) (numeric vector).
#' @return Leakage factor (m) (numeric vector).
#' @details The length of both input vectors should be the same.
.create_l <- function(kD, c) {
  sqrt(kD * c)
}

#' Internal function to calculate the quotient of the transmissivity and the leakage factor.
#'
#' @inheritParams .create_m
#' @param l Leakage factor (m) (numeric vector).
#' @return Quotient of the transmissivity and the leakage factor (m/d) (numeric vector).
#' @details The length of both input vectors should be the same.
.create_kDl <- function(kD, l) {
  kD / l
}

#' Internal function to calculate the matrix with the quotients of the radius of the areas and the leakage factors.
#'
#' @inheritParams .create_kDl
#' @param r Radius of areas (m) (numeric vector).
#' @return Matrix with the quotients of the radius of the areas and the leakage factors (-) (matrix).
#' @details The dimension of the matrix is (2, length(l)-1).
.create_rl <- function(r, l) {
  m <- length(l)
  i <- c(2:m)
  rbind(r / l[1:(m - 1)], r / l[i])
}

#' Internal function to create a zero two-row matrix.
#'
#' @param m The number of areas for which kD- and c-values are specified (-) (integer).
#' @return A zero two-row matrix. (matrix).
#' @details The dimension of the matrix is (2, length(l)-1).
.zero_rows <- function(m) {
  matrix(rep(0, 4 * m), nrow = 2)
}

#' Internal function to create two rows of the system matrix corresponding with area boundary j.
#'
#' @inheritParams .zero_rows
#' @param j area boundary index (-) (integer).
#' @param rl Matrix with the quotients of the radius of the areas and the leakage factors (-) (matrix).
#' @param kDl Quotient of the transmissivity and the leakage factor (m/d) (numeric vector).
#' @return Two rows of the system matrix corresponding with area boundary j (matrix).
#' @details The dimension of the matrix is (2, 2*m).
.create_rows <- function(j, m, rl, kDl) {
  res <- .zero_rows(m)
  i1 <- 1  + (j - 1) * 2
  i2 <- i1 + (m - 1) * 2 - 1
  res[1, i1:i2] <-
    c(besselI(rl[1, j], 0),
      besselK(rl[1, j], 0),-besselI(rl[2, j], 0),-besselK(rl[2, j], 0))
  res[2, i1:i2] <-
    c(
      kDl[j] * besselI(rl[1, j], 1),-kDl[j] * besselK(rl[1, j], 1),-kDl[j + 1] * besselI(rl[2, j], 1),
      kDl[j + 1] * besselK(rl[2, j], 1)
    )
  return(res)
}

#' Internal function to create the system matrix.
#'
#' @inheritParams .zero_rows
#' @inheritParams .create_kDl
#' @inheritParams .create_rows
#' @return System matrix (matrix).
#' @details The dimension of the matrix is (2m, 2*m).
.create_A_matrix <- function(m, l, kDl, rl) {
  res <- lapply(1:(m - 1), function(j)
    rbind(.create_rows(j, m, rl, kDl))) %>% do.call(rbind, .) %>% rbind(.zero_rows(m))
  res[2 * m - 1, 2] <- 1
  res[2 * m, 2 * m - 1] <- 1
  return(res)
}

#' Internal function to create the right hand vector b in the linear equations.
#'
#' @param h Polder heads in the different areas (m) (numeric vector).
#' @return right hand vector b in the linear equations.
#' @details The length of the vector h should be equal to the length of the vectors kD and c.
.create_b_vector <- function(h) {
  diff(h) %>% lapply(function(x)
    c(x, 0)) %>% unlist() %>% append(c(0, 0))
}

#' Internal function to determine the area number corresponding to a radius x.
#'
#' @inheritParams .create_rl
#' @param x Radius (m) (numeric).
#' @return Area number corresponding to a radius x (-) (integer).
.get_index <- function(x, r) {
  r <- c(0, r)
  stats::approx(
    r,
    y = 1:length(r),
    xout = x,
    method = "constant",
    rule = c(2:2)
  )$y
}

#' Initialise a list (y) with parameters used for the calculations of the groundwater heads and fluxes.
#'
#' @inheritParams .create_m
#' @inheritParams .create_l
#' @inheritParams .create_rl
#' @inheritParams .create_b_vector
#' @return List with parameters used for the calculations of the groundwater heads and fluxes.
#' @examples
#' kD <- c(1000, 1000, 1000)
#' c <- c(1000, 2000, 3000)
#' r <- c(1000, 2000)
#' h <- c(10, 9, 8)
#' y <- rd_init(kD, c, r, h)
#' @export
rd_init <- function(kD, c, r, h) {
  m <- .create_m(kD)
  l <- .create_l(kD, c)
  kDl <- .create_kDl(kD, l)
  rl <- .create_rl(r, l)
  A <- .create_A_matrix(m, l, kDl, rl)
  b <- .create_b_vector(h)
  coeff <- solve(A, b)
  list(
    m = m,
    l = l,
    kDl = kDl,
    rl = rl,
    c = c,
    r = r,
    h = h,
    coeff = coeff
  )
}

#' Internal function to create a list of the parameters corresponding to a distance x to be used for the calculations of the groundwater heads and fluxes (list).
#'
#' @inheritParams .get_index
#' @param y List with parameters used for the calculations of the groundwater heads and fluxes.
#' @return List of the parameters corresponding to a radius x to be used for the calculations of the groundwater heads and fluxes.
.get_coeff <- function(x, y) {
  i <- .get_index(x, y$r)
  list(
    Ai = y$coeff[1 + 2 * (i - 1)],
    Bi = y$coeff[2 + 2 * (i - 1)],
    rl = max(x,.001) / y$l[i],
    h = y$h[i],
    kDl = y$kDl[i],
    c = y$c[i]
  )
}

#' Calculate the head at radius x.
#'
#' @inheritParams .get_coeff
#' @return Head at radius x (m) (numeric).
#' @examples
#' kD <- c(1000, 1000, 1000)
#' c <- c(1000, 2000, 3000)
#' r <- c(1000, 2000)
#' h <- c(10, 9, 8)
#' y <- rd_init(kD, c, r, h)
#' x <- 500
#' rd_phi(x, y)
#' @export
rd_phi <- function(x, y) {
  y <- .get_coeff(x, y)
  y$Ai * besselI(y$rl, 0) + y$Bi * besselK(y$rl, 0) + y$h
}

#' Calculate the lateral discharge at radius x.
#'
#' @inheritParams .get_coeff
#' @return Lateral discharge (m3/d) (numeric).
#' @examples
#' kD <- c(1000, 1000, 1000)
#' c <- c(1000, 2000, 3000)
#' r <- c(1000, 2000)
#' h <- c(10, 9, 8)
#' y <- rd_init(kD, c, r, h)
#' x <- 500
#' rd_q(x, y)
#' @export
rd_q <- function(x, y) {
  y <- .get_coeff(x, y)
  - 2 * pi * y$kDl * x * (y$Ai * besselI(y$rl, 1) - y$Bi * besselK(y$rl, 1))
}

#' Internal function to determine if the radius x corresponds with a area boundary.
#'
#' @inheritParams .get_index
#' @inheritParams .create_rl
#' @return TRUE if radius x corresponds with a area boundary (Boolean).
.is_on_border <- function(x, r){
  which(r==x) %>% length() == 1
}


#' Calculate the seepage intensity at radius x.
#'
#' @inheritParams .get_coeff
#' @return Seepage intensity (m/d) (numeric).
#' @details Positive = downwards flow; negative = upwards flow.
#' @examples
#' kD <- c(1000, 1000, 1000)
#' c <- c(1000, 2000, 3000)
#' r <- c(1000, 2000)
#' h <- c(10, 9, 8)
#' y <- rd_init(kD, c, r, h)
#' x <- 500
#' rd_seep(x, y)
#' @export
rd_seep <- function(x, y) {
  if (!.is_on_border(x, y$r)) {
    phi <- rd_phi(x, y)
    y <- .get_coeff(x, y)
    (y$h - phi) / y$c
  } else {
    s1 <- rd_seep(x - .01, y)
    s2 <- rd_seep(x + .01, y)
    mean(c(s1,s2))
  }
}

#' Calculate the average seepage intensity between radius r1 and r2.
#'
#' @inheritParams .get_coeff
#' @param r1 Smallest radius (m) (numeric).
#' @param r2 Largest radius (m) (numeric).
#' @return Average seepage intensity (m/d) (numeric).
#' @details Positive = downwards flow; negative = upwards flow.
#' @examples
#' kD <- c(1000, 1000, 1000)
#' c <- c(1000, 2000, 3000)
#' r <- c(1000, 2000)
#' h <- c(10, 9, 8)
#' y <- rd_init(kD, c, r, h)
#' r1 <- 500
#' r2 <- 1500
#' av_seep(r1, r2, y)
#' @export
av_seep <- function(r1, r2, y) {
  Q1 <- rd_q(r1, y)
  Q2 <- rd_q(r2, y)
  (Q2 - Q1) / (pi * (r2^2 - r1^2))
}
