#' Test for the median in the one and two sample paired tests
#'
#' In the one sample case, test if the median of the random variable is equal to 0. In the paired two sample case, test if the median of the
#' difference between the two random variables is equal to 0.
#' @param x,y two continuous variables.
#' @param alternative indicates the alternative hypothesis and must be one of "two.sided", "greater" or "less".
#' @param paired a logical value. If it equals TRUE, you must provide values for \code{x} and \code{y}
#' and the paired test is implemented. If it equals FALSE, only \code{x} must be provided.
#' @details The null hypothesis for the one sample median test is: H0 Med(X)=0 where Med represents the median.
#' The alternative is specified by the \code{alternative} argument: "\code{greater}" means that Med(X)>0 and "\code{less}"
#'  means that Med(X)<0. The null hypothesis for the paired median test is: H0 Med(X-Y)=0. Both tests are asymptotically
#'  calibrated in the sense that the rejection probability under the null hypothesis is asymptotically equal to the level of the test. The
#'  test is based on the asymptotic law of the empirical median and uses a kernel estimator to estimate the density of \code{X} (in the one sample case)
#'  or of \code{X-Y} in the two sample case at 0.
#' @note The paired median test can be implemented by providing the variables \code{x} and \code{y} or by just providing
#' one vector equal to the difference between \code{x} and \code{y}.
#' @return Returns the result of the test with its corresponding p-value and the value of the test statistic.
#' @keywords test
#' @seealso \code{\link{cortest}}, \code{\link{indeptest}}, \code{\link{vartest}}, \code{\link{wilcoxtest}}.
#' @export
#' @examples
#' #Simulations
#' n=100
#' M=1000 #number of replications
#' res1=res2=rep(NA,M)
#' testone=function(n){
#' D=rchisq(n,df=4)-qchisq(df=4, p=0.5)
#' list(test1=mediantest(D)$p.value,test2=binom.test(sum(D>0),n)$p.value)
#' } #test2 is the sign test.
#' for (i in 1:M)
#' {
#' result=testone(n)
#' res1[i]=result$test1
#' res2[i]=result$test2
#' }
#' mean(res1<0.05) #0.048
#' mean(res2<0.05) # 0.04

mediantest <- function(x,y=NULL,alternative="two.sided",paired=FALSE) {UseMethod("mediantest")}
#' @export
mediantest.default=function(x,y=NULL,alternative="two.sided",paired=FALSE)
{
  if (paired==TRUE)
  {
    if (is.null(y)){ stop("'y' is missing for paired test")}
    if (is.null(x)){ stop("'x' is missing for paired test")}
    if (length(x)!=length(y)){stop("'x' and 'y' should have the same size for the paired median test")}
    #Perform the paired two sample test
    X <- x-y
    pairedTest<-TRUE
  }
  if (paired==FALSE)
  {
    if (is.null(x)) stop("'x' is missing")
    if (is.null(y)==FALSE) stop("there should be no 'y' for the one sample median test")
    X<-x
    pairedTest<-FALSE
  }
  n <- length(X)
  z <- stats::density(X, kernel="rectangular",n=2000)
  xab <- which.min(abs(z$x))
  medX=stats::median(X)
  Tn <- 2*z$y[xab]*sqrt(n)*medX
  if (alternative=="two.sided" | alternative=="t"){
    Pval <- 2*(1-pnorm(abs(Tn)))}
  if (alternative=="less"| alternative=="l"){
    Pval <- pnorm(Tn)}
  if (alternative=="greater"| alternative=="g"){
    Pval <- 1-pnorm(Tn)}
  result <- list(statistic=Tn, p.value=Pval,estimate=medX, alternative=alternative,pairedTest=pairedTest)
  class(result)<-"mediantest"
  return(result)
}

#' @export
print.mediantest <- function(x, ...)
{
  medval=x$estimate
  if (x$pairedTest==FALSE){
    cat("\nMedian test for the on sample case\n\n")
    cat(paste("Tn = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))
    if (x$alternative=="two.sided" | x$alternative=="t"){
      cat("alternative hypothesis: median (X) is not equal to zero\n")}
    if (x$alternative=="less" | x$alternative=="l"){
      cat("alternative hypothesis: median (X) is negative\n")}
    if (x$alternative=="greater" | x$alternative=="g"){
      cat("alternative hypothesis: median (X) is positive\n")}
    cat("median estimate:\n")
    print(medval)
  }
  if (x$pairedTest==TRUE){
    cat("\nMedian test for the two sample case\n\n")
    cat(paste("Tn = ", round(x$statistic,4), ", " , "p-value = ",round(x$p.value,4),"\n",sep= ""))
    if (x$alternative=="two.sided" | x$alternative=="t"){
      cat("alternative hypothesis: median (X-Y) is not equal to zero\n")}
    if (x$alternative=="less" | x$alternative=="l"){
      cat("alternative hypothesis: median (X-Y) is negative\n")}
    if (x$alternative=="greater" | x$alternative=="g"){
      cat("alternative hypothesis: median (X-Y) is positive\n")}
    cat("median estimate of the difference:\n")
    print(medval)
  }
}




