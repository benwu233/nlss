############### PNC data  ###################
library(gplots)
library(nlss)
library(gtools)
library(abind)
#################

"12 Dorsal attention"
"11 Ventral attention"
"3 Cingulo-opercular Task Control"
"8 Fronto-parietal Task Control"
"9 Salience"
"6 Memory retrieval"
"5 Default mode"
"10 Subcortical"
"13 Cerebellar"
"7 Visual"
"4 Auditory"
"1 Sensory/somatomotor Hand"
"2 Sensory/somatomotor Mouth"

source("/Users/Ben/Desktop/Rpackage/nlss_function.R")
group0 = read.csv("/Users/ben/Desktop/work2/data/PNC/roi_id_power.csv")

group = group0
order0 = c(12,11,3,8,9,6,5,10,13,7,4,1,2)
for(i in 1:13){
  group$network_id[group$network_id==order0[i]] = i+ 100
}
group$network_id[group$network_id!=-1] = group$network_id[group$network_id!=-1] - 100

group_1 = group$network_id
group_2 = group_1[order(group_1)][-c(1:28)]

FC_data_85 = readRDS(file = "/Users/ben/Desktop/work2/nlss/FC_data_new_power236_ordered_248s.rds")



set.seed(611)

res1 = NLSS(data=FC_data_85, q=8, kk = 1, group_node = group_2,
            total_iter = 3000, burn_in = 0, thin = 1, show_step=100)

meanS0 = res1$S[,,3000] -2
path0 = "/Users/ben/desktop/test/"
heatmap.net.new.nolg( meanS0,lim=c(-1,1),path = path0,community = group_2,color = bluered(3),filename = paste0("FC3030") )


S0 = res1$S[,,3000]
A0 = res1$A[,,3000]

simX = simnlss(A0,S0,3)

test = list()
test$simX = simX
test$A = A0
test$S = S0

saveRDS(test,"/Users/ben/desktop/test/test.rds")
test= readRDS("/Users/ben/desktop/test/test.rds")

res62 = NLSS(data=simX,init=res6, q=6, kk = 1, group_node = group_2,
            total_iter = 5000, burn_in = 0, thin = 1, show_step=100)

saveRDS(res,"/Users/ben/desktop/test/simres8.rds")
saveRDS(res2,"/Users/ben/desktop/test/simres8_2.rds")
saveRDS(res3,"/Users/ben/desktop/test/simres8_3.rds")

res = readRDS("/Users/ben/desktop/test/simres8.rds")
saveRDS(res6,"/Users/ben/desktop/test/simres6.rds")

res6 = readRDS("/Users/ben/desktop/test/simres6.rds")

init0 = list()
init0$S = res6$S[,,300]
init0$A = res6$A[,,300]
init0$beta = res6$beta[,,300]

init0 = list()
#init0$S = test$S
init0$A = test$A
#init0$beta = res3$beta[,,3000]

res7 = NLSS(data=test$simX,init=init0, q=8, kk = 1, group_node = group_2,
             total_iter = 3000, burn_in = 0, thin = 10, show_step=100)

loglik03 = log_lik_A(test$simX,res5,1,300,3)
loglik04 = log_lik_A(test$simX,res6,1,300,3)
loglik05 = log_lik_A(test$simX,res7,1,300,3)


plot(c(loglik03,loglik04,loglik05),ylim= c(-4465000,-4440000))
lines(rep(logliktrue,900),col="red" )
res_sum = NLSS_sum(res,2001,3000)

cor(t(res_sum$S),t(S0))

sum(res_sum$S[1,]!=S0[1,])

loglik0 = log_lik_A(test$simX,res,1,3000,3)
loglik00 = log_lik_A(test$simX,res2,1,3000,3)
loglik01 = log_lik_A(test$simX,res3,1,3000,3)
loglik02 = log_lik_A(test$simX,res4,1,3000,3)

plot(c(loglik0,loglik00,loglik01,loglik02),ylim= c(-4465000,-4440000))
lines(rep(logliktrue,12000),col="red" )

logliktrue = log_lik_A0(test$simX,test$S,test$A,3)




loglik1 = log_lik_A(test$simX,res6,1,3000,3)
loglik11 = log_lik_A(test$simX,res62,1,3000,3)
plot(loglik01)




meanS0 = init0$S-2
path0 = "/Users/ben/desktop/test/"
heatmap.net.new.nolg( meanS0,lim=c(-1,1),path = path0,community = group_2,color = bluered(3),filename = paste0("FC3232") )





path0 = paste0("/Users/ben/Desktop/work2/nlss/nlsstestpower236/convergence/")

A_list = readRDS(paste0(path0,"A_list_FC248.rds"))
log_list = readRDS(paste0(path0,"log_list_FC248.rds"))
beta_coef_list0 = readRDS(paste0(path0,"beta_coef_list_FC248.rds"))
S_list2 = readRDS(paste0(path0,"S_list2_FC248.rds"))

A_list = readRDS(paste0(path0,"A_list_FCSC.rds"))
log_list = readRDS(paste0(path0,"log_list_FCSC.rds"))
beta_coef_list0 = readRDS(paste0(path0,"beta_coef_list_FCSC.rds"))
S_list2 = readRDS(paste0(path0,"S_list2_FCSC.rds"))

meanS0_0 = list()
for(i in 1:5){
  betaco0 = beta_coef_list0[[i]]
  meanS095_0 = find_th(betaco0,0.95,2)
  meanS0_0[[i]] = meanS095_0 - 2
  #heatmap.net.new.nolg( meanS0_0[[i]],lim=c(-1,1),path = path1,community = group_2,color =  bluered(3), filename = paste0("FCold",i) )
}


meanS095_0 = find_th(beta_coef_list0[[1]],0.95,2) -2

rep0 = repro(meanS095_0,S_list2[[1]]-2,group_2)
L = length(group_2)
groupmat = matrix(0,L,L)

for(i in unique(group_2)){
  for(j in unique(group_2)){
    if(j>=i){
      groupmat[group_2==i,group_2==j] = i*100+j
      groupmat[group_2==j,group_2==i] =  i*100+j
    }
  }
}
group_vec = vec_mat(groupmat)

path2 = paste0("/Users/ben/Desktop/work2/nlss/nlsstestpower236/bs_match/")
beta_coef_list = readRDS(paste0(path2,"beta_coef_list_FC.rds"))
meanS0 = list()
for(i in 1:50){
  betaco0 = beta_coef_list[[i]]
  meanS095 = find_th(betaco0,0.95,1)
  meanS0[[i]] = meanS095 - 2
}


Stotal95=list()
Stotal95[[1]] = meanS095_0
for(i in 1:50){
  Stotal95[[i+1]] = meanS0[[i]]
}
S_match = match_rows(Stotal95)


relia00 = relia_rows(S_match,conn_only=FALSE)


S_list3 = list()
n=50
for(i in 1:5){
  Strace = S_list2[[i]] - 2
  Strace0 = array(0,c(8,27730,60))
  for(j in 1:60){
    print(paste(i,j))
    Stmp = Strace[,,((j-1)*n+1):(j*n) ]
    Strace0[,,j] = apply(Stmp,c(1,2),mean)
  }
  S_list3[[i]] = Strace0
}

group_unique = unique(group_vec)
L = length(group_unique)
net_list = list()
for(i in 1:5){
  net_trace = array(0,c(8,L,60))
  Strace0 = S_list3[[i]]
  for(l in 1:L){
    print(paste(i,l))
    net_trace[,l,] = apply(Strace0[,group_vec==group_unique[l],],c(1,3),mean)
  }
  net_list[[i]] = net_trace
}

plot(net_trace[1,64,])

path0 = "/Users/ben/desktop/test/"
for(i in 1:5){
  heatmap.net.new.nolg(meanS0_0[[i]],lim=c(-1,1),path = path0,community = group_2,color = bluered(3),filename = paste0("FC_test",i) )
}




