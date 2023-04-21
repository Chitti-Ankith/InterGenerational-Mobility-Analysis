#' Calculates 14 Intergenerational Mobility Indices including single stage and transition matrix indices
#'
#' Takes the SES data for Parents and their children, and prints out the value of 14 different mobility indices. Transition Matrix is also constructed by dividing the input data into several percentiles
#'
#' @param X Parents SES data (Preferably income)
#' @param Y Childs SES data (Preferably income)
#' @param grps Number of percentiles to divide the data into for calculating the transiton matrix
#' @param showtmat Determines whether to print the transition matrix or not
#' @export
#' @author Chitti Ankith Reddy,Bheemeshwar Reddy
#' @importFrom knitr kable
#' @importFrom utils install.packages
#' @importFrom stats cor cor.test ecdf lm quantile
#' @references Savegnago, Marco. (2016). Igmobil: A Command for Intergenerational Mobility Analysis in Stata. Stata Journal. 16. 386-402.

Mobility = function(X,Y,grps = 4,showtmat = FALSE){
  #Library "knitr" is required to print the data in the appropriate manner
  if(requireNamespace("knitr", quietly = TRUE)){
    print("knitr is loaded correctly")
  } else {
    print("trying to install knitr")
    install.packages("knitr")
    if(requireNamespace("knitr", quietly = TRUE)){
      print("knitr installed and loaded")
    } else {
      stop("could not install knitr")
    }
  }
  #X is the Income data for the parents
  #Y is the Income data for the children
  #grps is the required number of quantiles to divide into
  #showmat indicates whether to print the transition matrix or not
  if(length(X) != length(Y)){
    stop("Error : Data Column lengths don't match")
  }

  else{

    cat("\n\t\tSingle Stage Indices")
    s_no = c(1:10)
    #Names of the single stage indices
    indices = c("1/N * sum abs(X - Y)","1/N * sum (X - Y)^2","1/N * sum abs(ln X - ln Y)","1/N * sum abs(X/mu(X) - Y/mu(Y))",
                "1 - Pearson coef. (on logs)","1 - Spearman coef. (on logs)","1 - OLS(Y,X)","1 - OLS(ln Y,ln X)",
                "1/N * sum abs(CDF X - CDF Y)","1/N * sum (CDF X - CDF Y)^2")
    #Values of the respective single stage indices
    results = c(sum(abs(X - Y))/length(X),sum((X - Y)^2)/length(X),sum(abs(log(X) - log(Y)))/length(X),
                sum((X/mean(X) - Y/mean(Y))^2)/length(X), 1 - cor(log(X),log(Y)),1-cor.test(log(X),log(Y),method = "spearman",exact = FALSE)$estimate,
                1 - lm(Y ~ X)$coefficient[2],1 - lm(log(Y) ~ log(X))$coefficient[2],sum(abs(ecdf(Y)(Y) - ecdf(X)(X)))/length(X),
                sum((ecdf(X)(X) - ecdf(Y)(Y))^2)/length(X))
    #Creating the data frame for the single stage indices
    ssdf = data.frame(Sno = s_no,Index = indices,Value = results)
    print(kable(ssdf))

    #Dividing the data into the required number of quantiles based on the 'grps'input
    perc_data = data.frame(parent = X,child = Y)
    perc_data$pp = factor(cut(perc_data$parent, quantile(perc_data$parent,prob = seq(0,1,length = grps+1)), include.lowest = TRUE),labels = LETTERS[1:grps])
    perc_data$pp = as.numeric(perc_data$pp)
    perc_data$cp = factor(cut(perc_data$child, quantile(perc_data$child,prob = seq(0,1,length = grps+1)), include.lowest = TRUE),labels = LETTERS[1:grps])
    perc_data$cp = as.numeric(perc_data$cp)
    tmat = matrix(nrow = grps,ncol = grps)
    #Constructing the transition matrix based on conditional probability
    for(i in 1:grps)
    {
      for(j in 1:grps)
      {
        tmat[i,j] = mean(perc_data$cp[perc_data$pp == i] == j)
      }
    }

    #Only prints the Transition Matrix if the user asks to do so
    if(showtmat == TRUE)
    {
      cat("\n\tTransition Matrix\n")
      print(tmat)
    }
    #Calculating the eigen values of the transition matrix
    ev = eigen(tmat)
    k = ncol(tmat)
    cat("\n\tTransition Matrix Indices")
    #Names of the transition matrix indices
    tmatind = c("Shorrock/Prais index","Bartholomew index","1 - Second Largest Eigen Value","Determinant Index")
    s_no = c(11:14)
    bart = 0
    #bart is used for calculating the Bartholomew index
    for(i in 1 : nrow(tmat))
    {
      for(j in 1 : ncol(tmat))
      {
        bart = bart + tmat[i,j]*abs(i-j)
      }
    }
    #Results for the transition matrix indices
    tmatresults = c((k - sum(Re(ev$values)))/(k-1),bart/(k*(k-1)),1-abs(Re(ev$values[2])),1-prod(Re(ev$values)))
    #Creating the data frame for the transition matrix indices
    tmatdf = data.frame(Sno = s_no,Index = tmatind,Value = tmatresults)
    print(kable(tmatdf))
  }
}
