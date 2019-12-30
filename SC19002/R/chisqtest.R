#' @title  chi-square test
#' @description give two sample vectors,we do the chi-squre test with the null hypothesis that the two samples are connected, compute the list of chi-square test statistic,degree of freedom and p-value.
#' @importFrom stats pchisq complete.cases
#' @param x a numeric vector sample to be tested
#' @param y a numeric vector sample to be tested, the length of y should be equal to x's
#' @return a list of chi-square test statistic, degree of freedom and p-value.
#' @examples
#' \dontrun{
#' a <- c(1,3,4,5,4,7);
#' b <- 1:6;
#' chisqtest(a,b)
#' }
#' @export
chisqtest<- function(x, y){
  if (!is.numeric(x)) {
    stop("x must be numeric")}
  if (!is.numeric(y)) {
    stop("y must be numeric")}
  if (length(x) != length(y)) {
    stop("x and y must have the same length")}
  if (length(x) <= 1) {
    stop("length of x must be greater one")}
  if (any(c(x, y) < 0)) {
    stop("all entries of x and y must be greater or equal zero")}
  if (sum(complete.cases(x, y)) != length(x)) {
    stop("there must be no missing values in x and y")}
  if (any(is.null(c(x, y)))) {
    stop("entries of x and y must not be NULL")}
  m <- rbind(x, y)
  rs <- rowSums(m)
  cs <- colSums(m)
  n <- sum(m)
  me <- tcrossprod(rs, cs) / n
  stat = sum((m - me)^2 / me)
  df <- (length(rs) - 1) * (length(cs) - 1)
  p <- pchisq(stat, df = df, lower.tail = FALSE)
  return(list(statistic = stat, df = df, `p-value` = p))
}
