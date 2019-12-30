#' @title the sample skewness coefficient
#' @description This function is used to compute the skewness coefficient of a set of numeric data x
#' @param x the data which is a numeric vector
#' @return A numeric vector of length one
#' @examples
#' \dontrun{
#' x <- rnorm(100)
#' sk(x)
#' }
#' @export
sk <- function(x){
  xmean <- mean(x)
  m3 <- mean((x-xmean)^3)
  m2 <- mean((x-xmean)^2)
  m3/m2^1.5
}
