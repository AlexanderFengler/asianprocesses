X = as.data.table(matrix(rnorm(100,mean=0,sd=1),ncol=10))

X[,list(names(X)[1:3])]


X = data.table(X1 = c(1,1,1), X2 = c(2,2,2), V3 = c(1,2,3))
Y = data.table(Y1 = c(4,4,4), Y2 = c(5,5,5), V3 = c(1,2,3))

X[,sum(1),b=="V3",with=FALSE]


temp  = c_samples
setkeyv(temp,names(temp)[mu_s_len+1])
temp[,lapply(.SD,mean),by=key(temp),.SDcols = 1:2]
DT[,lapply(.SD,sum),by=x]
c_posterior_by_datapoint = c_samples[, list(p_c_1 = sum(names(c_posterior_by_datapoint)[i+1])/(sum(V1)+sum(V2)),p_c_2 = sum(V2)/(sum(V1)+sum(V2))),by = V3]
c_samples[,lapply(.SD,my_fun),by = V3,.SDcols = 1:2]


my_loop_1 = function(n){
  for (i in seq_len(n)){
  }
}

my_loop_2 = function(n){
  for (i in 1:n){
  }
}

microbenchmark::microbenchmark(my_loop_1(1000000),my_loop_2(1000000),times=100)


x = 1:4
y = 1:4

p = recordPlot()
plot.new()
p <- plot(x,y)
p

mu_samples_new = data.table(t(mu_samples),sample = 1:length(mu_samples[1,]))
c_samples_new = c_samples
setnames(c_samples_new,names(c_samples_new)[mu_s_len+2],'sample')

setkey(c_samples_new, sample)
setkey(mu_samples_new, sample)

c_samples_new[mu_samples_new]


cluster_order = apply(abs(matrix(rep(mu_s,each = length(mu_s)),nrow = length(mu_s), byrow = TRUE) - matrix(rep(apply(mu_samples,mean,MARGIN = 1),times = length(mu_s)), nrow = length(mu_s), byrow = TRUE)), which.min, MARGIN = 2)
c_samples_new = c_samples[cluster_order,]

dnorm(x_s %*% matrix(rep(1,mu_s_len*x_s_len),), mean = ((mu_0/sd_0^2) + (c_temp %*% x_s) / sd_s^2)  / ((1/sd_0^2) + (apply(c_temp, sum, MARGIN = 1)/sd_s^2)) %*% c_temp,
      sd = sqrt(1/((1/sd_0^2) + (apply(c_temp, sum, MARGIN = 1)/sd_s^2)) + sd_s^2),
      log = T)


