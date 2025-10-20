# Modelling-Claim-Severities-Using-the-Pareto-Lomax-Distribution-in-R
The analysis demonstrates a peaks-over-threshold (POT) approach for modeling heavy-tailed insurance claim severities using the Pareto Type II (Lomax) distribution — a standard tool in actuarial risk modeling and extreme value theory.
The claim data is simulated, select a high threshold using the mean excess function, fit a Lomax distribution via maximum likelihood and validate tail fit using log–log survival and Pareto-scale QQ diagnostics. Finally, key risk measures are computed which are; Value at Risk (VaR) and Tail Value at Risk (TVaR) — at the 99% confidence level.

> #Modelling Claim Severities using Pareto distribution
> #package for actuarial and risk modelling
> library(actuar) # Pareto/Lomax functions
> library(fitdistrplus)  # fitting univariate + diagnostics
> library(evir) # Extreme value theorem
> library(ggplot2) # Visualization 
> # Sample claims
> set.seed(1) # reproducibility
> # Type I (actuar::rpareto) distribution
> claims <- rpareto(3000, shape = 1.7, scale = 200) 
> 
> # Pick threshold using mean excess plot
> meplot(claims)  # choose u where curve looks ~linear in tail
> u <- quantile(claims, 0.80)  # 80th percentile of claims data
> 
> # Excesses over u (Lomax fit on excesses)
> excess <- claims[claims > u] - u
> 
> #  Fit Lomax (Pareto II) via MLE
> lomax_loglik <- function(par, y){
+   alpha <- par[1]; beta <- par[2]
+   if (alpha <= 0 || beta <= 0) return(Inf)
+   -sum(log(alpha/beta) - (alpha+1)*log1p(y/beta))
+ }
> #Evaluates the negative log-likelihood of the Lomax distribution
> #fit Lomax via optim
> fit <- optim(c(alpha=1.5, beta=median(excess)), lomax_loglik, y=excess, method="L-BFGS-B",
+              lower=c(1e-6, 1e-6))
> #maximum likelihood estimates of Lomax distribution’s parameters
> alpha_hat <- fit$par["alpha"]; beta_hat <- fit$par["beta"]
> #  Diagnostics: log-log survival plot
> library(ismev) # Extreme value modelling
> #Emphirical tail sample of exceedance over threshold
> emp_tail <- sort(excess, decreasing=TRUE) # Descending order
> #Empirical survival probabilities
> S_emp <- seq_along(emp_tail)/length(emp_tail)
> #Log-survival plot and Adjusted claim size
> plot(log(emp_tail+beta_hat), log(S_emp), main="Log-Log Tail", xlab="log(x+beta)", ylab="log S(x)")
> # QQ plot on Pareto scale
> z_emp <- sort(log1p(excess/beta_hat))
> # since log(1 + X/beta) ~ Exp(alpha)
> z_the <- qexp(ppoints(length(z_emp)), rate=alpha_hat) 
> qqplot(z_the, z_emp, xlab="Theoretical Exp Quantiles", ylab="Observed log(1+x/beta)")
> # straight Line passin through intercept 0 and slope of 1
> abline(0,1) 
> #  Risk metrics
> #Computed fitted VaR at probaility level p
> VaR <- function(p) beta_hat * ((1-p)^(-1/alpha_hat) - 1)
> #Computes TVaR for Pareto II distribution
> TVaR <- function(p) (alpha_hat*beta_hat)/(alpha_hat-1) * (1-p)^(-1/alpha_hat) - beta_hat
> VaR_99 <- VaR(0.99) # 99% VaR
> TVaR_99 <- TVaR(0.99) #99% TVaR
> #Creating the summary of results
> c(alpha_hat=alpha_hat, beta_hat=beta_hat, VaR_99=VaR_99, TVaR_99=TVaR_99)
alpha_hat.alpha   beta_hat.beta     VaR_99.beta   TVaR_99.alpha 
         1.6954        504.5015       7125.5091      18097.6274
>
> <img width="532" height="308" alt="image" src="https://github.com/user-attachments/assets/73ca6698-8849-48a6-8168-dcab0b684b6a" />
<img width="532" height="308" alt="image" src="https://github.com/user-attachments/assets/02467d5f-0bd2-49ac-9ca8-b9bce1c23862" />


> 
