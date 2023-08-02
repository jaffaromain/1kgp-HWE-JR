
# One-Sample Pearson Chi-Square Test
library(tidyverse)

Genotype_counts<-function(Naa,NAa,NAA){
  total <- Naa + NAa +  NAA
  #  genotype frequencies
  PAA = NAA/total 
  Paa = Naa/total
  PAa = NAa/total
  #observed allele frequencies
  Pa<-Paa  + (PAa/2) 
  PA<-1 - Pa
  Genotype_counts <- as.data.frame(cbind(PAA, PAa, Paa, Pa,PA,total))
  return(Genotype_counts)
}

Pearson_test <- function(Naa, NAa, NAA){
  genotype_counts <- Genotype_counts(Naa,NAa,NAA)
  #expected genotype frequencies
  genotype_counts <- genotype_counts %>% mutate(expected_AA = (PA^2) * total,
  expected_Aa = (2 * PA* Pa * total),expected_aa = (Pa^2) * total)
  
  n <- genotype_counts$total
  p <- genotype_counts$Pa
  # (Oaa - Eaa)^2/ Eaa
  E_aa = n*(p^2)
  E_Aa = n*(2*p*(1-p))
  E_AA <- n*((1-p)^2)
  stat = ((Naa - E_aa)^2/E_aa) + 
    ((NAa - E_Aa)^2/E_Aa) +
    ((NAA - E_AA)^2/E_AA)
  # get p-values
  vals <- c()
  for(i in 1:length(stat)){
    vals[i] <-  pchisq(stat[i], df =1, lower.tail = F)}
  genotype_counts$Pearson_p_val <- vals
  genotype_counts$T_Pearson <- stat
  return(data.frame(genotype_counts$T_Pearson, genotype_counts$Pearson_p_val, ref=p))
}

#########################################################################################################################
#Robust-allele test
Delta_Test <- function(Naa, NAa, NAA){
  genotype_counts <- Genotype_counts(Naa,NAa,NAA)
  # allele frequencies
  p=genotype_counts$Pa 
  n=genotype_counts$total 
  # get delta
  delta = (genotype_counts$Paa - p^2)
  delta_squared = delta^2 
  denominator = (1/n) * p^2 * (1-p)^2
  #get se
  se = sqrt(denominator)
  stat_delta <- delta_squared/denominator
  genotype_counts$T_delta <- stat_delta
  #iterate through all SNPs
  vals_delta <- c()
  for(i in 1:length(stat_delta)){
    vals_delta[i] <-  pchisq(stat_delta[i], df =1, lower.tail = F)}
  genotype_counts$Delta_p_val <- vals_delta
  genotype_counts$T_Delta <- stat_delta
  #store results
  return(data.frame(T_delta=genotype_counts$T_delta, Delta_p_val=genotype_counts$Delta_p_val,N=n, Paa=genotype_counts$Paa, 
                    PAa=genotype_counts$PAa, PAA=genotype_counts$PAA, PA=genotype_counts$PA ,Pa=genotype_counts$Pa ,delta, se))
}

####################  META-ANALYSIS SCRIPT  ################

#directional meta
meta_v1_inverse_variance<-function(results){
  #get delta and se
  beta_i=as.matrix(results%>%select(starts_with("Delta_")))
  se_i=as.matrix(results%>%select(starts_with("se")))
  # add weighting
  w_i= (1/se_i)^2
  se=sqrt(1/rowSums(w_i))
  beta=rowSums(beta_i*w_i)/rowSums(w_i)
  Z=beta/se
  meta_p=2*pnorm(abs(-Z),lower.tail = F)
  #  95%  confidence  interval
  L_CI=beta-se*1.96;U_CI=beta+se*1.96
  return(data.frame(SNP=results$SNP,Z=Z,`Meta P-value`=meta_p,L_CI,U_CI,beta,meta1_se=se))
}

#omnibus test
meta_v2<-function(results){
  # requires test statistic (delta/se^2)
  #get delta and se
  
  delta_i=as.matrix(results%>%select(starts_with("Delta_")))
  se_i=as.matrix(results%>%select(starts_with("se")))
  stat_delta <- (delta_i/se_i)^2
  test_stat=rowSums(stat_delta)# chisquare stat for each SNP
  #p-value
  p<- pchisq(test_stat, df=ncol(stat_delta), lower.tail = F)
  return(data.frame(SNP=results$SNP,Stat=test_stat,`Meta v2 P-value`=p))
}

# cochran's test
heterogen_test<-function(results){
  #get delta and se
  
  delta_i=as.matrix(results%>%select(starts_with("Delta_")))
  se_i=as.matrix(results%>%select(starts_with("se")))
  #test stat for SNP in each group
  stat_delta_z<- (delta_i/se_i)
  #overall test stat across all groups
  Q_stat=rowSums(stat_delta_z^2)-1/(rowSums((1/se_i)^2))*(rowSums(1/se_i*stat_delta_z))^2
  p <- pchisq(Q_stat, df=(ncol(stat_delta_z)-1), lower.tail = F)
  return(data.frame(SNP=results$SNP,`Q Het. P-value`=p))
}



