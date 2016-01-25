#' Metropolis version of the chinese restaurant process
#' @return
#' @export
#' @examples
fin_chin_c_mu_metro = function(n_data = 90, mu_s = c(-5,0,5), p_mu_s = c(0.6,0.2,0.2), sd_s = 0.5, sd_0 = 3, mu_0 = 0, alpha = 90/3, n_samples = 2000){

  library(gtools)
  library(data.table)
  library(hyperdirichlet)
  library(ggplot2)

  # GENERATE DATA --------------------------------------------------------
  # distribution x_i | c_i ~ N(mu_ci, 1)
  n_data = (n_data %/% (length(mu_s))) * length(mu_s)
  x_s = rep(0,rep(n_data)) # initialize vector for final samples
  temp = permute(rep(c(1:length(mu_s)),times=p_mu_s*length(x_s)))
  c_s = matrix(0,nrow = length(mu_s),ncol = length(x_s),byrow = 1) # we sue fixed categories

  for (i in 1:length(x_s)){
    c_s[temp[i],i] = 1
  }

  k = length(mu_s)# initialize dirichtlet parameter

  # mu_s = t(rnorm(n=3, sd = sd_0, mean = mu_0)) # generate three mu_1,..m_3 ~ N(mu_0,sd_0)
  # theta_s = rdirichlet(n = 1, alpha = rep(alpha,length(mu_s))) # generate thetas from dirichlet
  # c_s = rmultinom(length(x_s),size = 1,prob = theta_s) # generate a sample of class assignments

  mean_s = mu_s %*% c_s # generate a vector of means corresponding to class asssignments
  x_s = rnorm(length(x_s), mean = mean_s, sd = sd_s) # generate data sample x_s ~ N(mean_s, 1)
  # ----------------------------------------------------------------------

  # SAMPLING PROCESS -----------------------------------------------------
  # intialization of main variables
  x_s_len = length(x_s)
  mu_s_len = length(mu_s)

  c_temp = rmultinom(length(x_s),size = 1,prob = rep(1/length(mu_s),length(mu_s))) # Initialize random category sample
  P_c_I_c_x = c_temp # initialize P( c_i = k | c_-i, X )
  P_c_I_c_others = rep(0,length(mu_s)) # initialize P( c_i = k | c_-i )
  P_x_I_c = rep(0,length(mu_s)) # initialize P( X | c )

  mu_samples = matrix(rep(0,n_samples*length(mu_s)), ncol = n_samples, nrow = 3, byrow = TRUE)
  mu_samples[,1] = rnorm(n = mu_s_len,
                         mean = mu_0,
                         sd = sd_0)
  mu_prop = mu_samples[,1]

  c_samples = matrix(0,nrow = length(x_s)*n_samples,ncol = length(mu_s) + 2) # stores the final samples of category distributions
  c_samples[,length(mu_s) + 1] = rep(1:length(x_s),n_samples)
  c_samples[,length(mu_s) + 2] = rep(1:n_samples,each = length(x_s))
  c_samples[1:x_s_len,1:mu_s_len] = t(c_temp) # store current category sample

  loglik_post = rep(0,n_samples)
  category_cnt = apply(c_temp, sum, MARGIN = 1) + (alpha/k)

  loglik_post[1] = sum(dnorm(x_s,mean = mu_samples[,1] %*% c_temp, sd = sd_s,log = T)) +
    sum(dnorm(mu_samples[,1],mean = mu_0, sd = sd_0, log = T)) +
    (log(diri_norm(category_cnt)) - log(diri_norm(rep(alpha/k,mu_s_len))))

  x_seq = seq_len(x_s_len)
  cnt = 1
  new_sample = 0
  plot_update_ratio = n_samples %/% 20
  par(mfrow = c(1,1))

  while (cnt < n_samples){
    new_sample = 0
    mu_prop = mu_prop + rnorm(n = mu_s_len, mean = mu_0, sd = 0.2)

    for (i in x_seq){ # we go through each data point i
      temp_category_count = apply(c_temp[,-i], sum, MARGIN = 1) # sum of data points in each category excluding data point in consideration
      P_c_I_c_others = (temp_category_count + (alpha/x_s_len)) / (x_s_len - 1 + alpha) # fill in probability to belong to category i according to P(c_i = k | c_-i) = (m_i,k + a/K) / (N - 1 + a)
      P_x_I_c = dnorm(x_s[i], mean = mu_prop, sd = sd_s) # calculate p(X | c_i = k)
      P_c_I_c_x[,i] = P_c_I_c_others * P_x_I_c # get P(c_i = k | c_-i, X) as vector for each k
      c_temp[,i] = rmultinom(n = 1, size = 1, prob = P_c_I_c_x[,i]) # sample category for x_i (P_c_I_c_x is automatically normalize by rmultinom)
    }

    category_cnt = apply(c_temp, sum, MARGIN = 1) + (alpha/k)

    loglik_prop = sum(dnorm(x_s,mean = mu_prop %*% c_temp, sd = sd_s,log = T)) +
      sum(dnorm(mu_prop, mean = mu_0, sd = sd_0, log = T)) +
      (log(diri_norm(category_cnt)) - log(diri_norm(rep(alpha/k,mu_s_len))))
    acceptance_ratio = exp(loglik_prop - loglik_post[cnt])

    if (runif(1) < acceptance_ratio){
      cnt = cnt + 1
      loglik_post[cnt] = loglik_prop
      c_samples[((cnt-1)*x_s_len + 1):((cnt)*x_s_len),1:mu_s_len] = t(c_temp) # store current category sample
      mu_samples[,cnt] = mu_prop
      new_sample = 1
      print(cnt)
    } else {
      mu_prop = mu_samples[,cnt]
      c_temp = t(c_samples[((cnt-1)*x_s_len + 1):((cnt)*x_s_len),1:mu_s_len])
    }

    # Plot chain progress ----
    if (cnt %% plot_update_ratio == 0 && new_sample == 1){
      if (cnt < 2*plot_update_ratio){
        plot(1:cnt,loglik_post[1:cnt],
             type='l',
             ylim = c(min(loglik_post[(cnt %/% 2):cnt]) - abs(range(loglik_post[(cnt%/%2):cnt])[1] - range(loglik_post[(cnt%/%2):cnt])[2]),
                      max(loglik_post[(cnt%/%2):cnt]) + abs(range(loglik_post[(cnt%/%2):cnt])[1] - range(loglik_post[(cnt %/% 2):cnt])[2])),
             xlim = c(1,n_samples),
             ylab = 'Sample Log Likelihood',
             xlab = 'Sample Index')
        plot_cnt = cnt #+ 1
      } else {
        lines(plot_cnt:cnt,loglik_post[plot_cnt:cnt])
        plot_cnt = cnt #+ 1
      }

      print(plot_cnt)
      print(cnt)
    }
    # ----
  }
  # ----------------------------------------------------------------------

  # Reshuffle category names ---------------------------------------------
  c_samples_means = apply(mu_samples,mean,MARGIN = 1)
  cluster_order = c(apply(abs(matrix(rep(mu_s,each = length(mu_s)), nrow = length(mu_s), byrow = TRUE) -
                                matrix(rep(apply(mu_samples,mean,MARGIN = 1),times = length(mu_s)), nrow = length(mu_s), byrow = TRUE)), which.min, MARGIN = 1),
                    mu_s_len + 1, mu_s_len + 2)
  c_samples = c_samples[,cluster_order]
  # ----------------------------------------------------------------------


  # Sample proportions by category, by data point, from sampling process -
  category_sums = apply(c_samples[,1:mu_s_len],sum,MARGIN = 2)
  category_proportions_posterior = prop.table(category_sums)
  c_samples = as.data.table(c_samples)
  setkeyv(c_samples,names(c_samples)[mu_s_len+1])
  c_posterior_by_datapoint = c_samples[,lapply(.SD,mean),by=key(c_samples),.SDcols = 1:mu_s_len]
  # ----------------------------------------------------------------------

  # PLOTTING -------------------------------------------------------------
  # DATA
  par(mfrow = c(2,2))
  dat = data.table(data = x_s, category = as.vector(t(1:length(mu_s) %*% c_s)))

  plot(dat$data,rep(0,length(dat$data)),type='n',ylim=c(0,0.75),ylab = 'Density', xlab = 'Data Points',xlim=c(min(mu_s) - 4*sd_s, max(mu_s) + 4*sd_s), main = 'Data')
  for (i in unique(dat$category)){
    d = density(dat[category == i,data])
    polygon(d,col=i,border=i)
    rug(dat[category == i,data],col=i)
  }

  category_proportions_s = prop.table(tabulate(dat$category))
  barplot(category_proportions_s, names.arg = 1:length(category_proportions_s),
          ylab = 'Proportion', xlab = 'Category',
          col = adjustcolor(1:length(category_proportions_s),alpha.f = 0.3),border = NA,ylim = c(0,1),
          main = 'Data')

  # SAMPLING PROCESS
  plot(x_s,1:x_s_len,type='n',col = 1:mu_s_len %*% c_s,ylab='Index',xlab='Data Points',xlim=c(min(mu_s) - 4*sd_s, max(mu_s) + 4*sd_s))
  abline(h=1:length(x_s), col='grey20')
  points(x_s,1:x_s_len,pch='|',col = 1:mu_s_len %*% c_s)

  for (i in 1:mu_s_len){
    points(rep(mu_s[i],length(x_s)),1:x_s_len,pch=16,col = adjustcolor(i,alpha.f = 0.3),cex = t(c_posterior_by_datapoint[,i+1,with=FALSE])*1.5)
  }

  barplot(category_proportions_posterior, names.arg = 1:length(category_proportions_posterior),
          ylab = 'Proportion', xlab = 'Category',
          col = adjustcolor(1:length(category_proportions_posterior),alpha.f = 0.3),border = NA, ylim=c(0,1),
          main = 'Samples')

  par(mfrow = c(1,1))
  hist(mu_samples, breaks = n_samples*n_data / 500, main = 'Sampling Distribution of mu', xlab = 'Values', xlim=c(min(mu_s) - 4*sd_s, max(mu_s) + 4*sd_s))
  # ----------------------------------------------------------------------

  return(list(c_samples = c_samples,
              mu_samples = mu_samples,
              loglik_samples =  loglik_post,
              category_posteriors = c_posterior_by_datapoint,
              data = x_s,
              category_proportions_data = category_proportions_s))
}
