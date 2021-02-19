aft_bart_sp_HTE <- function(){
  # True median survival time
  true_mst_trt_1 <- rep(NA, dim(mydata$X)[1])
  for (i in 1:dim(mydata$X)[1]){
    fun_trt_1 <-
      function(x) {
        exp(-(1000/mydata$lambda1 * exp(mydata$LP1[i]) * x) ^ mydata$eta_1)-0.5
      }
    true_mst_trt_1[i] <- uniroot(fun_trt_1, interval = c(0, 100))$root
    # print(i)
  }
  
  true_mst_trt_2 <- rep(NA, dim(mydata$X)[1])
  for (i in 1:dim(mydata$X)[1]){
    fun_trt_2 <-
      function(x) {
        exp(-(1000/mydata$lambda2 * exp(mydata$LP2[i]) * x) ^ mydata$eta_2)-0.5
      }
    true_mst_trt_2[i] <- uniroot(fun_trt_2, interval = c(0, 100))$root
    # print(i)
  }
  
  
  # Create counterfactual data
  mydata_trt_1 <- mydata$X %>% as_tibble() %>%  
    filter(W == 1) %>% 
    select(-W) 
  mydata_trt_0 <- mydata$X %>% as_tibble() %>%  
    filter(W == 0) %>% 
    select(-W) 
  # Build ABART
  abart_mydata_trt_1 = abart(x.train = mydata_trt_1 %>% select(-Time,-Event) %>% as.matrix(), times = mydata_trt_1$Time, delta=mydata_trt_1$Event, nskip=100, ndpost=1000,x.test = mydata_trt_0 %>% select(-Time,-Event) %>% as.matrix())
  
  abart_mydata_trt_0 = abart(x.train = mydata_trt_0 %>% select(-Time,-Event) %>% as.matrix(), times = mydata_trt_0$Time, delta=mydata_trt_0$Event, nskip=100, ndpost=1000, x.test = mydata_trt_1 %>% select(-Time,-Event) %>% as.matrix())
  
  
  # Calculate the number of time points in ABRT output
  k_abart_mydata_trt_1 <- length(abart_mydata_trt_1$surv.train.mean)/dim(mydata_trt_1)[1]
  
  k_abart_mydata_trt_0 <- length(abart_mydata_trt_0$surv.train.mean)/dim(mydata_trt_0)[1]
  
  # Transform ABRT predicted survival vetor into subject-specific matrix
  abart_mydata_trt_1_trt_1_survival <- matrix(NA, nrow = dim(mydata_trt_1)[1], ncol = k_abart_mydata_trt_1)
  abart_mydata_trt_1_trt_0_survival <- matrix(NA, nrow = dim(mydata_trt_0)[1], ncol = k_abart_mydata_trt_1)
  abart_mydata_trt_0_trt_0_survival <- matrix(NA, nrow = dim(mydata_trt_0)[1], ncol = k_abart_mydata_trt_0)
  abart_mydata_trt_0_trt_1_survival <- matrix(NA, nrow = dim(mydata_trt_1)[1], ncol = k_abart_mydata_trt_0)
  
  for (i in 1:dim(mydata_trt_1)[1]){
    abart_mydata_trt_1_trt_1_survival[i,] <- abart_mydata_trt_1$surv.train.mean[(1+ k_abart_mydata_trt_1 * (i-1)):(k_abart_mydata_trt_1 + k_abart_mydata_trt_1 * (i-1))]
    abart_mydata_trt_0_trt_1_survival[i,] <- abart_mydata_trt_0$surv.test.mean[(1+ k_abart_mydata_trt_0 * (i-1)):(k_abart_mydata_trt_0 + k_abart_mydata_trt_0 * (i-1))]
    # print(i)
  }
  for (i in 1:dim(mydata_trt_0)[1]){
    abart_mydata_trt_1_trt_0_survival[i,] <- abart_mydata_trt_1$surv.test.mean[(1+ k_abart_mydata_trt_1 * (i-1)):(k_abart_mydata_trt_1 + k_abart_mydata_trt_1 * (i-1))]
    abart_mydata_trt_0_trt_0_survival[i,] <- abart_mydata_trt_0$surv.train.mean[(1+ k_abart_mydata_trt_0 * (i-1)):(k_abart_mydata_trt_0 + k_abart_mydata_trt_0 * (i-1))]
    # print(i)
  }
  
  abart_mydata_trt_1_survival <- abart_mydata_trt_1_trt_1_survival %>% 
    rbind(abart_mydata_trt_1_trt_0_survival)
  abart_mydata_trt_0_survival <- abart_mydata_trt_0_trt_1_survival %>% 
    rbind(abart_mydata_trt_0_trt_0_survival)
  
  # Calculate the time points in ABRT output 
  mydata_trt_1_time <- quantile(mydata_trt_1$Time, probs = c(1:k_abart_mydata_trt_1)/k_abart_mydata_trt_1)
  mydata_trt_0_time <- quantile(mydata_trt_0$Time, probs = c(1:k_abart_mydata_trt_0)/k_abart_mydata_trt_0)
  # Get the predicted median surviavl time from ABART model
  predict_abart_mydata_median_survival_trt_1 <- rep(NA, dim(mydata$X)[1])
  predict_abart_mydata_median_survival_trt_0 <- rep(NA, dim(mydata$X)[1])
  for (j in 1:dim(mydata$X)[1]){
    predict_abart_mydata_median_survival_trt_1_index <- which(abs(abart_mydata_trt_1_survival[j,]-0.5) == min(abs(abart_mydata_trt_1_survival[j,]-0.5)))
    
    predict_abart_mydata_median_survival_trt_1[j] <- mydata_trt_1_time[predict_abart_mydata_median_survival_trt_1_index]
    
    predict_abart_mydata_median_survival_trt_0_index <- which(abs(abart_mydata_trt_0_survival[j,]-0.5) == min(abs(abart_mydata_trt_0_survival[j,]-0.5)))
    
    predict_abart_mydata_median_survival_trt_0[j] <- mydata_trt_0_time[predict_abart_mydata_median_survival_trt_0_index]
    
  }
  
  # Test
  bias_ate <- tibble(p1_category = 1:50,
                     true_ate = rnorm(50, 14.1,0.5),
                     bias_ate = rnorm(50,0.28,0.05))
  
  # Calculate the bias
  true_ate <- tibble(true_mst_trt_2 = true_mst_trt_2, 
                     true_mst_trt_1 = true_mst_trt_1,
                     p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(true_ate = true_mst_trt_2 -  true_mst_trt_1) %>% 
    summarise(true_ate = mean(true_ate)) %>% 
    ungroup()
  
  abart_ate <- tibble(predict_mst_trt_2 = predict_abart_mydata_median_survival_trt_1, 
                      predict_mst_trt_1 = predict_abart_mydata_median_survival_trt_0,
                      p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(ate = predict_mst_trt_2 -  predict_mst_trt_1) %>%
    # pull(ate) %>% mean
    summarise(ate = mean(ate)) 
  
  Bias_ate <- true_ate %>% 
    inner_join(abart_ate) %>% 
    mutate(bias = true_ate - ate)
  
  return(bias_ate)
}


