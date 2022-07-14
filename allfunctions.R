#Rcpp::sourceCpp("GUniFrac.cpp")



"%w/o%" <- function(x, y) x[!x %in% y] #--  x without y

### microbial distances BC
UniFrac_BC = function(OTUtab = OTUtab) {return(as.matrix(stats::dist(OTUtab, method = 'manhattan'))/2)}

### Hamming distances
UniFrac_Hamming = function(OTUtab = OTUtab,weight =1) {
  OTUtab_pre = OTUtab
  OTUtab_pre[OTUtab>0] = 1 # change abundance to presence information
  weight_OTUtab_pre = t(t(OTUtab_pre)*weight)
  return(as.matrix(stats::dist(weight_OTUtab_pre,method = "manhattan")))
}

###case-control matrix
case_control_mat = function(response) {
  # response: a vector of binary outcomes (0s and 1s)
  # Construct a matrix to show the case/control status between responses: 1 (same) -1 (different)
  return(1-2*as.matrix(stats::dist(response, method = 'manhattan')))
}

#######weight calculation for different rho###########
weight.cal<-function(rho,first_term,n_ker){
  new_w = rep(0,n_ker)
  if(rho==0){
    tempmax = which.max(first_term)
    new_w[tempmax]=1
  } else{
    new_w = exp((first_term-max(first_term))  / rho) # the new weights for kernels
    new_w = new_w/sum(new_w) # scale the kernels
  }
  return(new_w)  
}





##########boosted microbiome distances UniGrac########
weightedUni<- function (otu.tab, tree, cum, alpha = c(0, 0.5, 1),weight) 
{
  alpha = matrix(alpha, 1, length(alpha))
  #alpha2 = matrix(c(0, 0.5, 1), 1, 3)
  #GuniF.cum = cum #cumcouting计算，一直要用，可避免重复计算
  #cum = GuniF.cum$cum
  br.len <- tree$edge.length*weight
  br.len = as.matrix(br.len)
  
  GUniF = GUniFracCpp2(cum, br.len, alpha)
  output = list(GUniF = GUniF, alpha = alpha)
  
  output
}



cumcounting<- function (otu.tab, tree) 
{
  if (!is.rooted(tree)) 
    stop("Rooted phylogenetic tree required!")
  otu.tab <- as.matrix(otu.tab)
  row.sum <- rowSums(otu.tab)
  otu.tab <- otu.tab/row.sum
  n <- nrow(otu.tab)
  # Construct the returning array
  if (is.null(rownames(otu.tab))) {
    rownames(otu.tab) <- paste("comm", 1:n, sep="_")
  }
  
  # Check OTU name consistency
  if (sum(!(colnames(otu.tab) %in% tree$tip.label)) != 0) {
    stop("The OTU table contains unknown OTUs! OTU names
					in the OTU table and the tree should match!" )
  }
  
  # Get the subtree if tree contains more OTUs
  absent <- tree$tip.label[!(tree$tip.label %in% colnames(otu.tab))]
  if (length(absent) != 0) {
    tree <- drop.tip(tree, absent)
    warning("The tree has more OTU than the OTU table!")
  }
  tip.label <- tree$tip.label
  otu.tab <- otu.tab[, tip.label]
  edge <- tree$edge
  br.len <- tree$edge.length
  cum <-  cumall (otu.tab,edge,br.len)
  cum.ct <- round(t(t(cum) * row.sum))
  return(list(cum.counts=cum.ct,cum.abun=cum))
}

###########OTU weights calculation######

weight_chisq <- function(datama,ind1,ind0,n_train){
  cum_pre = datama
  cum_pre[datama>0] = 1 # change abundance to presence information
  cum_pre_table_11 = rowSums(cum_pre[,ind1])
  cum_pre_table_10 = rowSums(cum_pre[,ind0])
  cum_pre_table_01 = length(ind1)-cum_pre_table_11
  cum_pre_table_00 = length(ind0)-cum_pre_table_10
  cum_pre_table = cbind(cum_pre_table_11,cum_pre_table_10,cum_pre_table_01,cum_pre_table_00)
  weight_a_c = apply(cum_pre_table,1,function(x) mychisq(x,n_train))
  #weight_a_c[is.na(weight_a_c)] = 1
  return(weight_a_c)
}

mychisq <- function(x,n_train)
{
  if (sum(x[1:2])==0){
    p=1
  } else if (sum(x[1:2])==n_train) {
    p = 1
  } else {
    p = chisq.test(matrix(x,nrow = 2,byrow = TRUE))$p.value
  }
  return(p)
}

mytt <- function(x,ind1,ind0){
  tmp = unique(x)
  tmp2 = length(tmp)
  if(tmp2==1){
    return(1)
  }else if(isTRUE(all.equal(tmp,rep(1,tmp2))) ) {
    return(1)
  }else{
    obj<-try(t.test(x[ind1], y = x[ind0]), silent=TRUE)
    if (is(obj, "try-error")) return(NA) else  return(obj$p.value)
  }
  
}

#######cross-validation for rho#########
rho_cv <- function(Kmat_ls, CC, y, rho_ls=c(0,100,1000,10000), kfold=5){
  n_sam = length(y)
  n_ker = length(Kmat_ls) # number of kernels
  ## Split data for cross-validatoin
  foldid=sample(rep(seq(kfold), length=n_sam))
  
  output = matrix(0,nrow = kfold, ncol = length(rho_ls))
  
  for (i in 1:kfold) {
    test_ind_2 = which(foldid==i)
    tran_ind_2 = which(foldid!=i)
    train_Kmat_ls=list()
    for (ll in 1:n_ker) { 
      train_Kmat_ls[[ll]]=Kmat_ls[[ll]][tran_ind_2,tran_ind_2]#######training kernel matrix
    }
    
    train_CC=CC[tran_ind_2,tran_ind_2]
    
    first_term = vector() # 1st term in optimization
    for (ll in 1:n_ker) { 
      first_term[ll] = sum(train_Kmat_ls[[ll]] * train_CC)
    }
    ########
    #rho_ls=c(0,0.5,1,2,5,10)
    weight_ls<-sapply(rho_ls, function(y) weight.cal(y,first_term,n_ker))

    auc_ls<-apply(weight_ls,2, function(x) {tmp= pre_cal(x,n_sam,n_ker,Kmat_ls,y[tran_ind_2],tran_ind_2,test_ind_2);pROC::roc(y[test_ind_2], tmp, levels = c(0,1), direction = '<')$auc})
    output[i,] = auc_ls
  }
  tmpid = which.max(colMeans(output))
  return(rho_ls[tmpid])
  
}

##########searching max value of rho#######
rhomax = function( first_term, rho = 100, n_ker ) {
  
  thisRho = rho
  tol = 1e-4
  u = n_ker
  logU = log(u)
  
  # compute weight
  thisW = weight.cal(thisRho,first_term,n_ker)
  H = sum(-thisW*log(thisW))
  # evaluate whether the perplexity is within tolerance
  Hdiff = H - logU
  tries = 0
  while (abs(Hdiff) > tol && tries < 30) {
    #if not, increase or decrease precision
    thisRho = 2*thisRho
    # compute the new values
    thisW = weight.cal(thisRho,first_term,n_ker)
    H = sum(-thisW*log(thisW))
    # evaluate whether the perplexity is within tolerance
    Hdiff = H - logU
    tries = tries + 1
  }
  
  return(thisRho)
}





#######predicton outcome given weight#######
pre_cal<-function(new_w,n_sam,n_ker,Kmat_ls,Y,tran_ind,test_ind){
  ##############prediction###############
  
  w_s_Kmat = matrix(0, n_sam, n_sam) # weighted average of scaled kernel matrices 
  for (ll in 1:n_ker) {w_s_Kmat = w_s_Kmat + new_w[ll] * Kmat_ls[[ll]]}
  simi_all_domain=w_s_Kmat
  
  case_tran_ind =tran_ind[ which(Y==1)]
  ctrl_tran_ind = tran_ind[ which(Y==0)]
  n_train=length(tran_ind)
  Y_tran_pred_simi_t = rep(NA, n_train)
  
  for (i in 1:n_train) {
    idx = tran_ind[i]
    x = simi_all_domain[idx,]
    x = x[-idx]
    caseIndex = x[na.omit(case_tran_ind%w/o%idx)]
    ctrlIndex = x[na.omit(ctrl_tran_ind%w/o%idx)]
    Y_tran_pred_simi_t[i] = t.test(caseIndex, ctrlIndex)$statistic
  }
  
  simi_logst = glm(Y~Y_tran_pred_simi_t, family = binomial())
  
  
  n_valid=n_sam-n_train
  Y_test_pred_simi_t = rep(NA, n_valid)
  
  # for each testing patient
  for (i in 1:n_valid) {
    x = simi_all_domain[test_ind[i],]
    caseIndex = x[case_tran_ind]
    ctrlIndex = x[ctrl_tran_ind]
    Y_test_pred_simi_t[i] = t.test(caseIndex, ctrlIndex)$statistic
  }
  
  Y_test_pred_simi = predict(simi_logst, newdata = tibble(Y_tran_pred_simi_t=Y_test_pred_simi_t), type = "response")
  
  return(Y_test_pred_simi)
}


