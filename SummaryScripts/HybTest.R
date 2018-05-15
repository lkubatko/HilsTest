#sample1 is a matrix of results read in
#e.g.: sample1 = matrix(scan("results.txt"),ncol=15,byrow=T)
#invar = 1 returns test for first invariant
#invar = 2 returns test for second invariant
#invar = 3 returns test for sum of invariants

hybtest = function(sample1,invar) {

  n = sum(sample1[1,])  
  p9 = sample1[,9]/n
  p7 = sample1[,7]/n
  p4 = sample1[,4]/n
  p3 = sample1[,3]/n
  p2 = sample1[,2]/n
  p6 = sample1[,6]/n
  
  
  ## values for first pair of invariants
  obs_invp1.a = sample1[,9]-sample1[,7]
  obs_invp2.a = sample1[,4]-sample1[,7]
  obs_var_invp1.a = n*p9*(1-p9) + n*p7*(1-p7) + 2*n*p9*p7
  obs_var_invp2.a = n*p4*(1-p4) + n*p7*(1-p7) + 2*n*p4*p7
  obs_cov_invp1_invp2.a = -1*n*p9*p4 + n*p9*p7 + n*p7*p4 + n*p7*(1-p7) 
  ratio.a = obs_invp2.a/obs_invp1.a
  h_invp1.true = 0
  h_invp2.true = 1
  GH.ts.a = (obs_invp1.a)*(ratio.a - h_invp1.true/h_invp2.true)/sqrt(obs_var_invp1.a*(ratio.a^2)-2*obs_cov_invp1_invp2.a*ratio.a+obs_var_invp2.a)
  
  ## values for second pair of invariants
  obs_invp1.b = sample1[,3]-sample1[,2]
  obs_invp2.b = sample1[,3]-sample1[,6]
  obs_var_invp1.b = n*p3*(1-p3) + n*p2*(1-p2) + 2*n*p3*p2
  obs_var_invp2.b = n*p3*(1-p3) + n*p6*(1-p6) + 2*n*p3*p6
  obs_cov_invp1_invp2.b = n*p3*(1-p3) + n*p3*p6 + n*p3*p2 - n*p2*p6 
  ratio.b = obs_invp2.b/obs_invp1.b
  h_invp1.true = 0
  h_invp2.true = 1
  GH.ts.b = (obs_invp1.b)*(ratio.b - h_invp1.true/h_invp2.true)/sqrt(obs_var_invp1.b*(ratio.b^2)-2*obs_cov_invp1_invp2.b*ratio.b+obs_var_invp2.b)
  
  ## values for sum of the pairs of invariants
  obs_invp1.c = obs_invp1.a - obs_invp1.b
  obs_invp2.c = obs_invp2.a - obs_invp2.b
  obs_var_invp1.c = n*p9*(1-p9) + n*p7*(1-p7) + n*p3*(1-p3) + n*p2*(1-p2) +
                    2*n*p9*p7 + 2*n*p9*p3 - 2*n*p9*p2 - 2*n*p7*p3 + 2*n*p7*p2 + 2*n*p3*p2
  obs_var_invp2.c = n*p4*(1-p4) + n*p7*(1-p7) + n*p3*(1-p3) + n*p6*(1-p6) +
                    2*n*p4*p7 + 2*n*p4*p3 - 2*n*p4*p6 - 2*n*p7*p3 + 2*n*p7*p6 + 2*n*p3*p6
  obs_cov_invp1_invp2.c = -1*n*p9*p4 + n*p9*p7 + n*p9*p3 - n*p9*p6
                          +  n*p7*p4 + n*p7*(1-p7) - n*p7*p3 + n*p7*p6
                          +  n*p3*p4 - n*p7*p3 + n*p3*(1-p3) + n*p3*p6
                          -  n*p2*p4 + n*p2*p7 + n*p2*p3 - n*p2*p6
  ratio.c = obs_invp2.c/obs_invp1.c
  h_invp1.true = 0
  h_invp2.true = 1
  GH.ts.c = (obs_invp1.c)*(ratio.c - h_invp1.true/h_invp2.true)/sqrt(obs_var_invp1.c*(ratio.c^2)-2*obs_cov_invp1_invp2.c*ratio.c+obs_var_invp2.c)
  
  if (invar==1) return(GH.ts.a)
  else if (invar==2) return(GH.ts.b)
  else if (invar==3) return(GH.ts.c)
  }