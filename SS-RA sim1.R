
library(MASS)
library(susieR)
seed = 1245
n = 50000
p.snp.X = 25
p.snp.Y = 0
p.snp.U = 5
n0 = 500
rho = 0.5
# xi.range = 0.2
# theta.XtoY = 0.1
# theta.YtoX = 0

summary_sim_corr_adj <- function(seed, rho, xi.range, theta.XtoY, theta.YtoX = 0){
  set.seed(seed)  
  n = n # sample size
  n0 = n0
  p.snp.X = p.snp.X # number of SNPs in g_X
  p.snp.Y = p.snp.Y # number of SNPs in g_Y
  p.snp.U = p.snp.U # number of SNPs in g_B
  xi.range = xi.range # generate non-zero xi's from uniform distribution, implying correlated pleiotropy
  theta.XtoY = theta.XtoY # causal effect from X to Y
  theta.YtoX = theta.YtoX # causal effect from Y to X
  p.all = p.snp.X + p.snp.Y + p.snp.U
  
  alpha.effect = runif(p.snp.X,0.2,0.3)*(rbinom(p.snp.X,1,0.5)*2-1)
  
  beta.effect = runif(p.snp.Y,0.2,0.3)*(rbinom(p.snp.Y,1,0.5)*2-1)
  
  ### generate effects of invalid IVs
  gamma.effect = runif(p.snp.U,0.2,0.3)*(rbinom(p.snp.U,1,0.5)*2-1)
  eta.effect = runif(p.snp.U,0.2,0.3)*(rbinom(p.snp.U,1,0.5)*2-1)
  xi.effect = runif(p.snp.U,-xi.range,xi.range)
  
  Z_uncorrelated <- matrix(rbinom((2*n+n0)*p.all, 2, MAF), ncol = p.all, byrow = TRUE)
  p.correlated <- floor(p.all / 2) 
  cor_matrix <- diag(1, p.all)
  cor_matrix[1:p.correlated, 1:p.correlated] <- rho
  diag(cor_matrix) <- 1  
  L <- chol(cor_matrix)
  Z <- Z_uncorrelated %*% L
  U = Z%*%c(rep(0,p.snp.X),rep(0,p.snp.Y),xi.effect) + rnorm(2*n+n0,0,sqrt(2))
  error_X = rnorm(2*n+n0,0,1)
  error_Y = rnorm(2*n+n0,0,1)
  
  X = 
    (Z%*%c(alpha.effect,theta.YtoX*beta.effect,gamma.effect + theta.YtoX*eta.effect) + 
       (1+theta.YtoX)*U + error_X + theta.YtoX*error_Y)/(1-theta.XtoY*theta.YtoX)
  
  Y = 
    (Z%*%c(theta.XtoY*alpha.effect,beta.effect,theta.XtoY*gamma.effect + eta.effect) + 
       (1+theta.XtoY)*U + theta.XtoY*error_X + error_Y)/(1-theta.XtoY*theta.YtoX)
  
  ### Generate two independent samples (Z1,X1,Y1) and (Z2,X2,Y2)
  Z1 = Z[1:n,]
  X1 = as.matrix(X[1:n])
  Y1 = as.matrix(Y[1:n])
  
  Z2 = Z[(n+1):(2*n),]
  X2 = as.matrix(X[(n+1):(2*n)])
  Y2 = as.matrix(Y[(n+1):(2*n)])
  
  Z0 = Z[(2*n+1):(2*n + n0),]
  Z0 = scale(Z0,scale = F)
  
  Z1 = scale(Z1,scale = F)
  X1 = scale(X1,scale = F)
  Y1 = scale(Y1,scale = F)
  
  ### Get summary statistics for X from marginal linear regression
  z1Tz1 = t(Z1)%*%Z1
  b_X = as.numeric(1/diag(z1Tz1)*(t(Z1)%*%X1))
  rep_X1 = X1[,rep(1,p.all)]
  se_X = 
    sqrt(colSums((rep_X1 - Z1%*%diag(b_X))^2)/(n-2)/diag(z1Tz1))
  
  betaZX <- b_X
  se_betaZX <- se_X 
  
  ### Get summary statistics for Y from marginal linear regression
  Z2 = scale(Z2,scale = F)
  X2 = scale(X2,scale = F)
  Y2 = scale(Y2,scale = F)
  
  z2Tz2 = t(Z2)%*%Z2
  b_Y = as.numeric(1/diag(z2Tz2)*(t(Z2)%*%Y2))
  rep_Y2 = Y2[,rep(1,p.all)]
  se_Y = 
    sqrt(colSums((rep_Y2 - Z2%*%diag(b_Y))^2)/(n-2)/diag(z2Tz2))
  
  betaZY <- b_Y
  se_betaZY <- se_Y
  
  FirstStage = lm(X1~0+Z1)
  gamma_hat = FirstStage$coef
  gamma_hat[is.na(gamma_hat)] <- 0
  sigma2_hat = sigma(FirstStage)^2
  
  corr <- t(Z2)%*%Z2/n
  corr0 <- t(Z0)%*%Z0/n0
  
  a1 = estimate_s_rss(n = n2, z = betaZY/se_betaZY, corr0, r_tol = 1e-10, method = "null-mle")
  corr0.adj = (1-a1)*corr0 + (a1)*diag(dim(corr0)[1])
  
  Ix0 = sum_MR_SuSiE_sel(corr, betaZX, betaZY, se_betaZY, gamma_hat, n1, n2)
  Ix0
  result.est = stage2(length(Ix0), corr[Ix0, Ix0], betaZX[Ix0], betaZY[Ix0], 
                      se_betaZY[Ix0], n1 = n, n2 = n, gamma_hat[Ix0], sigma2_hat)
  p.est = result.est$Pvalue
  b.est = result.est$Estimate
  se.est = result.est$SE
  p.est
  b.est
  se.est
  
  Ix0 = sum_MR_SuSiE_sel(corr0.adj, betaZX, betaZY, se_betaZY, gamma_hat, n1, n2)
  Ix0
  result.est = stage2(length(Ix0), corr0.adj[Ix0, Ix0], betaZX[Ix0], betaZY[Ix0], 
                      se_betaZY[Ix0], n1 = n, n2 = n, gamma_hat[Ix0], sigma2_hat)
  p.est2 = result.est$Pvalue
  b.est2 = result.est$Estimate
  se.est2 = result.est$SE
  
  Ix0 = sum_MR_SuSiE_sel(corr0, betaZX, betaZY, se_betaZY, gamma_hat, n1, n2)
  Ix0
  result.est = stage2(length(Ix0), corr0[Ix0, Ix0], betaZX[Ix0], betaZY[Ix0], 
                      se_betaZY[Ix0], n1 = n, n2 = n, gamma_hat[Ix0], sigma2_hat)
  p.est3 = result.est$Pvalue
  b.est3 = result.est$Estimate
  se.est3 = result.est$SE
  
  Ix.oracle = 1:p.snp.X
  Ix.oracle
  result.oracle = stage2(length(Ix.oracle), corr[Ix.oracle, Ix.oracle], betaZX[Ix.oracle], betaZY[Ix.oracle], 
                         se_betaZY[Ix.oracle], n1 = n, n2 = n, gamma_hat[Ix.oracle], sigma2_hat)
  p.oracle = result.oracle$Pvalue
  b.oracle = result.oracle$Estimate
  se.oracle = result.oracle$SE
  
  result.oracle = stage2(length(Ix.oracle), corr0.adj[Ix.oracle, Ix.oracle], betaZX[Ix.oracle], betaZY[Ix.oracle], 
                         se_betaZY[Ix.oracle], n1 = n, n2 = n, gamma_hat[Ix.oracle], sigma2_hat)
  p.oracle2 = result.oracle$Pvalue
  b.oracle2 = result.oracle$Estimate
  se.oracle2 = result.oracle$SE
  
  result.oracle = stage2(length(Ix.oracle), corr0[Ix.oracle, Ix.oracle], betaZX[Ix.oracle], betaZY[Ix.oracle], 
                         se_betaZY[Ix.oracle], n1 = n, n2 = n, gamma_hat[Ix.oracle], sigma2_hat)
  p.oracle3 = result.oracle$Pvalue
  b.oracle3 = result.oracle$Estimate
  se.oracle3 = result.oracle$SE
 
  ## other methods
  ZTZ = corr
  ZTY = matrix(diag(ZTZ) * betaZY, ncol = 1);
  YTY = NULL
  for(SNP in 1:dim(ZTZ)[1]){
    YTY[SNP] = (n-1)*ZTZ[SNP, SNP]*(se_betaZY^2)[SNP] + ZTY[SNP]*betaZY[SNP]
  }
  YTY = median(YTY)
  
  # IWAS
  beta <- matrix(ginv(t(betaZX) %*% ZTZ %*% betaZX) %*% t(betaZX) %*% (ZTY), ncol = 1)
  sigma <- as.numeric(unlist((YTY - t(ZTY)%*%betaZX%*%ginv(t(betaZX)%*%ZTZ%*%betaZX)%*%t(betaZX)%*%ZTY)/(n - 1)))
  se <- matrix(sqrt(diag(ginv(t(betaZX) %*% ZTZ %*% betaZX) * c(sigma))), ncol = 1)
  Z <- matrix(beta/se, ncol = 1)
  p_xy_iwas <- matrix(2*pnorm(-abs(Z), 0, 1), ncol = 1)  
  
  # TWAS
  zGWAS = betaZY/se_betaZY
  z_uvIWAS = betaZX%*%zGWAS/sqrt(betaZX%*%ZTZ%*%betaZX)
  p_xy_twas <- as.numeric(2*pnorm(-abs(z_uvIWAS), 0, 1))
  
  return(list(
    pxy = c(p.est, p.est2, p.est3, p.oracle, p.oracle2, p.oracle3, p_xy_iwas, p_xy_twas)
  ))
}  




edge1 = matrix(data = NA, nrow = 200, ncol = 8)
edge2 = matrix(data = NA, nrow = 200, ncol = 8)
edge3 = matrix(data = NA, nrow = 200, ncol = 8)
edge4 = matrix(data = NA, nrow = 200, ncol = 8)
edge5 = matrix(data = NA, nrow = 200, ncol = 8)
edge6 = matrix(data = NA, nrow = 200, ncol = 8)
edge7 = matrix(data = NA, nrow = 200, ncol = 8)

for (i in 1:200){
  sim1 = summary_sim_corr_adj(seed = i, rho = 0.8, xi.range = 0.2, theta.XtoY = -0.05, theta.YtoX = 0.0)
  sim2 = summary_sim_corr_adj(seed = i, rho = 0.8, xi.range = 0.2, theta.XtoY = -0.03, theta.YtoX = 0.0)
  sim3 = summary_sim_corr_adj(seed = i, rho = 0.8, xi.range = 0.2, theta.XtoY = -0.01, theta.YtoX = 0.0)
  sim4 = summary_sim_corr_adj(seed = i, rho = 0.8, xi.range = 0.2, theta.XtoY = 0.0, theta.YtoX = 0.0)
  sim5 = summary_sim_corr_adj(seed = i, rho = 0.8, xi.range = 0.2, theta.XtoY = 0.01, theta.YtoX = 0.0)
  sim6 = summary_sim_corr_adj(seed = i, rho = 0.8, xi.range = 0.2, theta.XtoY = 0.03, theta.YtoX = 0.0)
  sim7 = summary_sim_corr_adj(seed = i, rho = 0.8, xi.range = 0.2, theta.XtoY = 0.05, theta.YtoX = 0.0)
  edge1[i,] = sim1$pxy
  edge2[i,] = sim2$pxy
  edge3[i,] = sim3$pxy
  edge4[i,] = sim4$pxy
  edge5[i,] = sim5$pxy
  edge6[i,] = sim6$pxy
  edge7[i,] = sim7$pxy
}

edge1 = na.omit(edge1)
edge2 = na.omit(edge2)
edge3 = na.omit(edge3)
edge4 = na.omit(edge4)
edge5 = na.omit(edge5)
edge6 = na.omit(edge6)
edge7 = na.omit(edge7)

edge1

pw1_xy = c(length(which(edge1[,1]<0.05))/dim(edge1)[1], length(which(edge2[,1]<0.05))/dim(edge2)[1],
           length(which(edge3[,1]<0.05))/dim(edge3)[1], length(which(edge4[,1]<0.05))/dim(edge4)[1],
           length(which(edge5[,1]<0.05))/dim(edge5)[1], length(which(edge6[,1]<0.05))/dim(edge6)[1],
           length(which(edge7[,1]<0.05))/dim(edge7)[1])

pw2_xy = c(length(which(edge1[,2]<0.05))/dim(edge1)[1], length(which(edge2[,2]<0.05))/dim(edge2)[1],
           length(which(edge3[,2]<0.05))/dim(edge3)[1], length(which(edge4[,2]<0.05))/dim(edge4)[1],
           length(which(edge5[,2]<0.05))/dim(edge5)[1], length(which(edge6[,2]<0.05))/dim(edge6)[1],
           length(which(edge7[,2]<0.05))/dim(edge7)[1])

pw3_xy = c(length(which(edge1[,3]<0.05))/dim(edge1)[1], length(which(edge2[,3]<0.05))/dim(edge2)[1],
           length(which(edge3[,3]<0.05))/dim(edge3)[1], length(which(edge4[,3]<0.05))/dim(edge4)[1],
           length(which(edge5[,3]<0.05))/dim(edge5)[1], length(which(edge6[,3]<0.05))/dim(edge6)[1],
           length(which(edge7[,3]<0.05))/dim(edge7)[1])

pw4_xy = c(length(which(edge1[,4]<0.05))/dim(edge1)[1], length(which(edge2[,4]<0.05))/dim(edge2)[1],
           length(which(edge3[,4]<0.05))/dim(edge3)[1], length(which(edge4[,4]<0.05))/dim(edge4)[1],
           length(which(edge5[,4]<0.05))/dim(edge5)[1], length(which(edge6[,4]<0.05))/dim(edge6)[1],
           length(which(edge7[,4]<0.05))/dim(edge7)[1])

pw5_xy = c(length(which(edge1[,5]<0.05))/dim(edge1)[1], length(which(edge2[,5]<0.05))/dim(edge2)[1],
           length(which(edge3[,5]<0.05))/dim(edge3)[1], length(which(edge4[,5]<0.05))/dim(edge4)[1],
           length(which(edge5[,5]<0.05))/dim(edge5)[1], length(which(edge6[,5]<0.05))/dim(edge6)[1],
           length(which(edge7[,5]<0.05))/dim(edge7)[1])

pw6_xy = c(length(which(edge1[,6]<0.05))/dim(edge1)[1], length(which(edge2[,6]<0.05))/dim(edge2)[1],
           length(which(edge3[,6]<0.05))/dim(edge3)[1], length(which(edge4[,6]<0.05))/dim(edge4)[1],
           length(which(edge5[,6]<0.05))/dim(edge5)[1], length(which(edge6[,6]<0.05))/dim(edge6)[1],
           length(which(edge7[,6]<0.05))/dim(edge7)[1])

pw7_xy = c(length(which(edge1[,7]<0.05))/dim(edge1)[1], length(which(edge2[,7]<0.05))/dim(edge2)[1],
           length(which(edge3[,7]<0.05))/dim(edge3)[1], length(which(edge4[,7]<0.05))/dim(edge4)[1],
           length(which(edge5[,7]<0.05))/dim(edge5)[1], length(which(edge6[,7]<0.05))/dim(edge6)[1],
           length(which(edge7[,7]<0.05))/dim(edge7)[1])

pw8_xy = c(length(which(edge1[,8]<0.05))/dim(edge1)[1], length(which(edge2[,8]<0.05))/dim(edge2)[1],
           length(which(edge3[,8]<0.05))/dim(edge3)[1], length(which(edge4[,8]<0.05))/dim(edge4)[1],
           length(which(edge5[,8]<0.05))/dim(edge5)[1], length(which(edge6[,8]<0.05))/dim(edge6)[1],
           length(which(edge7[,8]<0.05))/dim(edge7)[1])


beta = c(-0.05, -0.03, -0.01, 0, 0.01, 0.03, 0.05)

df1 = data.frame(methods=rep(c("MR-SuSiE-SS", "MR-SuSiE-SS-RA", "MR-SuSiE-SS-R", "Oracle-SS", "Oracle-SS-RA", "Oracle-SS-R", 'IWAS', 'TWAS'), each = 7),
                 beta, c(pw1_xy, pw2_xy, pw3_xy, pw4_xy, pw5_xy, pw6_xy, pw7_xy, pw8_xy))
colnames(df1) = c('methods', 'samplesize', 'values')

df1 = data.frame(methods=rep(c("MR-SuSiE-SS", "MR-SuSiE-SS-RA", "MR-SuSiE-SS-R"), each = 7),
                 beta, c(pw1_xy, pw2_xy, pw3_xy))
colnames(df1) = c('methods', 'samplesize', 'values')
df2 = data.frame(methods=rep(c("Oracle-SS", "Oracle-SS-RA", "Oracle-SS-R"), each = 7),
                 beta, c(pw4_xy, pw5_xy, pw6_xy))
colnames(df2) = c('methods', 'samplesize', 'values')

plot1 = ggplot(df1, aes(x=samplesize, y=values, color=methods)) + 
  geom_line(aes(linetype=methods), size = 1.5) +
  geom_point(aes(shape=methods), size = 3, stroke = 1.5) +
  scale_linetype_manual(values = c(6,7,8)) +
  scale_shape_manual(values=c(6,7,8)) + 
  scale_color_brewer(palette = "Set2") +
  scale_x_continuous(breaks = beta) + 
  ylab("Empirical power / type-I error") +
  xlab('theta') + 
  ylim(0,1) + 
  theme(text = element_text(size = 15))

plot2 = ggplot(df2, aes(x=samplesize, y=values, color=methods)) + 
  geom_line(aes(linetype=methods), size = 1.5) +
  geom_point(aes(shape=methods), size = 3, stroke = 1.5) +
  scale_linetype_manual(values = c(6,7,8)) +
  scale_shape_manual(values=c(6,7,8)) + 
  scale_color_brewer(palette = "Set2") +
  scale_x_continuous(breaks = beta) + 
  xlab("theta") + ylab("Empirical power / type-I error") +
  ylim(0,1) + 
  theme(text = element_text(size = 15))

grid.arrange(plot1, plot2, ncol = 2)


