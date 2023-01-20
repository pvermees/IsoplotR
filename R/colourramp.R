delta <- 6 / 29
delta2 <- delta * delta
delta3 <- delta2 * delta

curve <- function(x) {
  if (delta3 < x) x ^ (1 / 3)
  else x / (3 * delta2) + 4 / 29
}

anticurve <- function(y) {
  if (delta < y) y ^ 3
  else (y - 4 / 29) * (3 * delta2)
}

rgb2xyz <- matrix(c(
  0.49 / 95, 0.0018, 0,
  0.31 / 95, 0.0081, 0.01 / 109,
  0.2 / 95, 0.0001, 0.99 / 109
), ncol = 3)
xyz2rgb <- solve(rgb2xyz)

# R, g and b from 0 to 1
rgb2ceilab <- function(rgb) {
  v <- rgb2xyz %*% rgb
  curvex <- curve(v[1])
  curvey <- curve(v[2])
  curvez <- curve(v[3])
  return(c(
    116 * curvey - 16,
    500 * (curvex - curvey),
    200 * (curvey - curvez)
  ))
}

ceilab2rgb <- function(lab) {
  l <- lab[1]
  a <- lab[2]
  b <- lab[3]
  curvey <- (l + 16) / 116
  curved <- c(a / 500 + curvey, curvey, curvey - b / 200)
  xyz <- vapply(curved, anticurve, 0)
  xyz2rgb %*% xyz
}

clamp01 <- function(x) pmin(1, pmax(0, x))

rgbv <- function(v) {
  c <- clamp01(v)
  rgb(c[1], c[2], c[3])
}

colour_ramp <- function(
  colour_list,
  collapse_red_green = FALSE,
  collapse_blue_yellow = FALSE
) {
  collapse <- c(1,
    if (collapse_red_green) 0.001 else 1,
    if (collapse_blue_yellow) 0.001 else 1
  )
  last_colour <- colour_list[[length(colour_list)]]
  colour_count <- length(colour_list)
  rgb_matrix <- col2rgb(colour_list) / 255
  lab_list <- lapply(1:colour_count, function(i) {
    rgb2ceilab(rgb_matrix[, i]) * collapse
  })
  lab_diffs <- mapply(
    function(u, v) (u - v),
    lab_list[2:colour_count],
    lab_list[1:colour_count - 1]
  )
  dists <- unlist(lapply(
    1:(colour_count - 1),
    function(i) sqrt(sum(lab_diffs[, i]^2))
  ))
  dist_order <- order(dists, decreasing = TRUE)
  total_dist <- sum(dists)
  function(n) {
    max_dist <- total_dist / (n - 2)
    min_counts <- floor(dists / max_dist)
    total <- sum(min_counts)
    extras <- n - total - 1
    which_sides <- dist_order <= extras
    counts <- min_counts + which_sides
    # which colours are we interpolating for each result colour?
    i2coli <- rep(1:(colour_count - 1), times = counts)
    # how many colours are on each interpolation?
    i2count <- rep(counts, times = counts)
    isstart <- c(TRUE, i2coli[1:(n - 2)] != i2coli[2:(n - 1)])
    # how far are we interpolating for each result colour?
    i2sidei <- (1:(n - 1)) - cummax((1:(n - 1)) * isstart)
    results <- mapply(
      function(coli, sidei, count) {
        lab_list[[coli]] + sidei * lab_diffs[, coli] / count
      },
      i2coli,
      i2sidei,
      i2count,
      SIMPLIFY = FALSE
    )
    rgb_results <- lapply(results, function(v) {
      ceilab2rgb(v / collapse)
    })
    c(vapply(rgb_results, rgbv, ""), last_colour)
  }
}

colour_iridesce <- colour_ramp(
  list("black", "dark magenta", "dark blue", "green", "yellow"),
  collapse_red_green = TRUE
)

colour_salmon <- colour_ramp(
  list("black", "red", "pink"),
  collapse_red_green = TRUE
)

colour_rainbow_plus <- colour_ramp(
  list("black", "red", "yellow", "green", "blue", "magenta")
)

colour_red_blue_black <- colour_ramp(
  list("red", "magenta", "blue", "black"),
  collapse_red_green = TRUE
)
