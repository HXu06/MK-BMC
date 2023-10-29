########main function to do prediction#####
######
#subOTUtable: whole OTU table for both training and testing data
#subtree: phylogenetic tree
#y: binary outcome for training data
#Lcov: kernel list for covariates for both training and testing data
#tran_ind: index of training data in subOTUtable
#test_ind: index of testin data in subOTUtable
#rho_ls: given grid of rho_ls, default is NA

mklpre<-function(subOTUtable, subtree, y, Lcov,tran_ind, test_ind, rho_ls= NA, grid = 10, kfold=5){
  n_sam = nrow(subOTUtable)
  
  ind1 = which(y==1)
  ind0 = which(y==0)
  
  n_train=length(tran_ind)
  n_test=length(test_ind)
  
  data_train=subOTUtable[tran_ind,]
  data_test=subOTUtable[test_ind,]
  
  
  data_train_abun = data_train/rowSums(data_train) # convert OTU counts to abundance
  data_test_abun = data_test/rowSums(data_test) # convert OTU counts to abundance
  
  dataset_abun = subOTUtable/rowSums(subOTUtable) # convert OTU counts to abundance
  
  ##cum_for all
  cum_all = cumcounting(subOTUtable,subtree)
  
  
  ntip = NCOL(dataset_abun)
  tmp = subtree$edge[,2]
  tmp2 = order(tmp)[1:ntip]
  
  
  #########weighted distance##############
  weight_a_c = weight_chisq(cum_all$cum.counts[,tran_ind],ind1=ind1,ind0 = ind0,n_train)
  weight_l_c = weight_a_c[tmp2]
  weight_a_t = apply(cum_all$cum.abun[,tran_ind], 1,function(x) mytt(x,ind1 = ind1, ind0 = ind0))
  weight_l_t =weight_a_t[tmp2]
  
  
  weight_a_c = -log10(weight_a_c)
  weight_a_c = weight_a_c/sum(weight_a_c)
  weight_l_c = -log10(weight_l_c)
  weight_l_c = weight_l_c/sum(weight_l_c)
  weight_a_t = -log10(weight_a_t)
  weight_a_t = weight_a_t/sum(weight_a_t)
  weight_l_t = -log10(weight_l_t)
  weight_l_t = weight_l_t/sum(weight_l_t)
  
  weight_dataset_abun = t(t(dataset_abun)*weight_l_t)
  weight_D_BC = UniFrac_BC(weight_dataset_abun) # Bray-Curtis distance
  weight_D_U = weightedUni(otu.tab = subOTUtable, subtree,cum =  cum_all$cum.abun, alpha = 1,weight_a_c)$GUniF[,,2] # Unweighted UniFrac distance
  weight_D_W = weightedUni(otu.tab = subOTUtable, subtree, cum =  cum_all$cum.abun,alpha = 1,weight_a_t)$GUniF[,,1] # Weighted UniFrac distance
  weight_D_H = UniFrac_Hamming(subOTUtable,weight_l_c)  #Hamming distance
  
  
  CC=case_control_mat(y)
  diag(CC)=0
  ###########Kernel matrix list#####
  distmat_ls = list(weight_D_BC,weight_D_U,weight_D_W,weight_D_H ) # a list of distance matrices
  Kmat_ls= distmat_ls
  distmat_ls_2 = list()
  for (ll in 1:length(Kmat_ls)) {
    distmat_ls_2[[ll]] = distmat_ls[[ll]][tran_ind,tran_ind]
    sigma = mean(distmat_ls_2[[ll]][upper.tri(distmat_ls_2[[ll]])]) # mean of distances
    Kmat_ls[[ll]] = exp(- distmat_ls[[ll]]^2 / 2/ sigma^2 ) # tranform distances to kernels
    
  }
  names(distmat_ls) = names(Kmat_ls) = c("W-BC","W-Unweighted","W-weighted","W-Hamming")  
  
  Kmat_ls = c(Kmat_ls,list(Lcov))
  #########################
  
  n_ker = length(Kmat_ls) # number of kernels

  train_Kmat_ls=list()
  for (ll in 1:n_ker) { 
    train_Kmat_ls[[ll]]=Kmat_ls[[ll]][tran_ind,tran_ind]#######training kernel matrix
  }
  
  train_CC=CC[tran_ind,tran_ind]
  
  first_term = rep(0,n_ker) # 1st term in optimization
  for (ll in 1:n_ker) { 
    first_term[ll] = sum(train_Kmat_ls[[ll]] * train_CC)
  }
  ########
  if(is.na(rho_ls)) {
    maxrho = rhomax( first_term, rho = 100, n_ker )
    rho_ls2 = seq(from=0, to=log10(maxrho), length.out = grid)
    rho_ls = c(0,10^rho_ls2)
  }
  
  
  rho.cv = rho_cv(train_Kmat_ls, train_CC, y, rho_ls, kfold)
  
  weight_ls<-sapply(rho.cv, function(y) weight.cal(y,first_term,n_ker))
  pre_ls<-apply(weight_ls,2, function(x) pre_cal(x,n_sam,n_ker,Kmat_ls,y,tran_ind,test_ind))
  output <- list(weight = weight_ls, outcome_pr = pre_ls, rho = rho.cv)
  
  return(output)
}
