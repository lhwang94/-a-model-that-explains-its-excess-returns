rm(list=ls())
# This function estimates the regression coefficients in multivariate
# linear regressions. It also computes the respective standard errors.

# Obtaining the coefficients estimates
OLS <- function(X,y) {
  
  beta=solve(t(X) %*% X) %*% t(X) %*% y 
  # Obtaining the fitted values from the regression
  yhat= X%*% beta         
  
  # Obtaining the number of observations and the number of regressors 
  # (including the the constant) employed in the analysis
  nobs=dim(as.matrix(X))[1]
  nreg=dim(as.matrix(X))[2]
  
  # Compute the vector of residuals
  res=y-yhat                           
  
  # Computing estimates of the error variance
  shat=as.numeric((t(res) %*% res)/(nobs-nreg))          
  
  # Computing the variance co-variance matrix of the beta estimates
  varcovar= shat * solve(t(X) %*% X)                
  
  # Computing simple t-statistic for the null of beta equal to zero
  tstat=beta/sqrt(diag(varcovar))
  
  # White's Standard Errors
  sandwich_var_covar <- solve(t(X) %*% X) %*%
    t(X) %*% diag(as.numeric(res^2)) %*% X %*%
    solve(t(X) %*% X)
  tsta_with_Wh = beta/sqrt(diag(sandwich_var_covar))
  
  # Computing the residual sum of squares
  RSS=t(res) %*% res             
  
  # Computing the total sum of squares
  TSS=t((y-mean(y))) %*% (y-mean(y)) 
  
  # Computing the R-squared
  Rq=1-(RSS/TSS)      
  
  # Computing the adjusted R-squared
  Rqad=1- ((nobs-1)/(nobs-nreg))*(1-Rq); 
  
  S=list( beta = beta, yhat = yhat, res = res, varcovarBeta = varcovar,
          varRes = shat, tstat = tstat, tsta_with_Wh = tsta_with_Wh, Rq = Rq, Rqadj = Rqad)
} 

# Load the library we need to load and save excel files. 
# Install if it does not load for you
library(readxl)
library(xlsx)
library(optimbase)
library(MASS)
library(data.table)

# Defining parameters of the code
# Determining the significance for the confidence interval
sig <- 0.05

# Loading monthly data for portfolio returns
data<-data.frame(read_excel("E:/UMD/BUFN640/factor_monthly.xlsx"))
class(data)
str(data)

# Isolating the non_dble returns under consideration
non_duble_ret <- data[,3]

# Isolating the risk-free rate vector
rf<-data[,9]

# Constructing excess returns for non_dble
mkt_exc <- data[,4]

# Constructing excess returns for the portfolios under consideration
non_duble_exc_ret<- non_duble_ret-rf

# Isolating the FF Factors
FF_factors <- data[,5:8]

#isolating the macroeconomic and other financial variables
macro_and_financial_var <- data[,10:21]

# Obtaining the number of observations for the non_dble returns
n <- dim(FF_factors)[1]

# Constructing the matrix of regressors
X <- as.matrix(cbind(rep(1,n),mkt_exc,FF_factors,macro_and_financial_var))
colnames(X)[1]<- "constant"

# Conducting the OLS analysis on the non_dble 
OLS_estimates <- OLS(X,non_duble_exc_ret)

# Obtaining the number of regressors
nreg <- dim(X)[2]
plot(OLS_estimates$res,type="l", col="blue", lwd=3)

# White's Test for Heteroskedasticity
sqr_res <- OLS_estimates$res^2

plot(sqr_res, type="l")

#backward elimination
FitAll <- lm(non_duble_exc_ret ~ Mkt.RF + SMB + HML + RMW + CMA + Divident + 
               d_p + PMI + M2_real + CPI_all + Emp_nondbles + Personal_Income + Housing + Retail_sales + 
               Ex_rate_avg + PPI + IP..cons.nondble, data = data)
summary(FitAll)
step(FitAll, direction = "backward")
lmbe=lm(formula= non_duble_exc_ret ~ Mkt.RF + SMB + RMW + Divident + 
          d_p + PMI + M2_real + CPI_all + Emp_nondbles + Personal_Income + 
          Housing + Ex_rate_avg + PPI, data = data)
summary(lmbe)

#test zero-mean shocks
exp_res <- (1/n)*sum(lmbe$residuals)
print(paste('expectation_shocks = ', 
            round(exp_res, digits=2)))

#test uncorrelated shocks
## Computing the Durbin-Watson test
library(lmtest)
dwtest(formula = non_duble_exc_ret ~ Mkt.RF + SMB + RMW + Divident + 
         d_p + PMI + M2_real + CPI_all + Emp_nondbles + Personal_Income + 
         Housing + Ex_rate_avg + PPI,data=data)

##test uncorrlationness between shocks by Computing Breusch-Godfrey Test with maximum lag equal to 3.
bgtest(formula = non_duble_exc_ret ~ Mkt.RF + SMB + RMW + Divident + 
         d_p + PMI + M2_real + CPI_all + Emp_nondbles + Personal_Income + 
         Housing + Ex_rate_avg + PPI, order = 3, data = data)

##test uncorrlationness between shocks by Computing Box-Pierce Test with 3-th order sample
Box.test(x = lmbe$residuals, lag = 3)

#test normality
##Jarque-Bera test of normality

xi_hat=(1/n)*sum((lmbe$residuals-mean(lmbe$residuals))^3)/
  ((1/n)*sum((lmbe$residuals-mean(lmbe$residuals))^2))^(3/2)


k_hat=(1/n)*sum((lmbe$residuals-mean(lmbe$residuals))^4)/
  ((1/n)*sum((lmbe$residuals-mean(lmbe$residuals))^2))^(2)

JarqueBera_test=n*((xi_hat^2)/6+((k_hat-3)^2)/24)
pval_JarqueBera_test=1-pchisq(JarqueBera_test,2)

print(paste('Jarque-Bera statistic = ', 
            round(JarqueBera_test, digits=2)))
print(paste('p-value of Jarque-Bera statistic = ', 
            round(pval_JarqueBera_test, digits=2)))

# Plotting the histogram and the fitted normal
fit <- fitdistr(lmbe$residuals, "normal")
para <- fit$estimate
hist(lmbe$residuals,breaks = 40, freq = F , 
     xlim=c(-0.5,0.5), ylab = 'Probability')
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)

#shift y to t+1
yx=shift(non_duble_ret, n=1, fill=NA, type="lead")
xyp=data.frame(cbind(X,yx))
#regressor selection for predicting model
formula= as.formula(paste("yx~", paste(colnames(xyp[,2:18]), collapse= "+")))
FitAll <- lm(formula=formula, data = xyp)
step(FitAll, direaction = "backward")
lmbep=lm(formula= yx ~ CMA + d_p + PMI + M2_real + Housing, data = xyp)
summary(lmbep)
pred=predict(FitAll,se.fit = TRUE)
pred$fit[length(pred$fit)]

#in December 2011, the predicted result is 0.009646399 , so BUY

library(lmtest)
#DW test and normality test
dwtest(formula= yx ~ CMA + d_p + PMI + M2_real + Housing,data=xyp)


xi_hat=(1/n)*sum((lmbe$residuals-mean(lmbe$residuals))^3)/
  ((1/n)*sum((lmbe$residuals-mean(lmbe$residuals))^2))^(3/2)


k_hat=(1/n)*sum((lmbe$residuals-mean(lmbe$residuals))^4)/
  ((1/n)*sum((lmbe$residuals-mean(lmbe$residuals))^2))^(2)

JarqueBera_test=n*((xi_hat^2)/6+((k_hat-3)^2)/24)
pval_JarqueBera_test=1-pchisq(JarqueBera_test,2)

print(paste('Jarque-Bera statistic = ', 
            round(JarqueBera_test, digits=2)))
print(paste('p-value of Jarque-Bera statistic = ', 
            round(pval_JarqueBera_test, digits=2)))

# Plotting the histogram and the fitted normal
fit <- fitdistr(lmbep$residuals, "normal")
para <- fit$estimate
hist(lmbep$residuals,breaks = 40, freq = F , 
     xlim=c(-0.5,0.5), ylab = 'Probability')
curve(dnorm(x, para[1], para[2]), col = 2, add = TRUE)