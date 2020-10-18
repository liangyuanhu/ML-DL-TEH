sl_lib_g <- c( "SL.randomForest","SL.gam") #choose your own esemble algorithm here 
sl_lib_censor <- c( "SL.randomForest","SL.gam")
sl_lib_failure <- c( "SL.randomForest","SL.gam")

TSHEE_HTE <- function(){
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
  
  mydata_df <- mydata$X %>% as.data.frame
  
  # Fit Superlearner
  sl_fit <- initial_sl_fit(
    T_tilde = mydata_df$Time,
    Delta = mydata_df$Event,
    A =  mydata_df$W,
    W = mydata_df[,1:10],
    t_max = max(mydata_df$Time),
    sl_treatment = sl_lib_g,
    sl_censoring = sl_lib_censor,
    sl_failure = sl_lib_failure
  )
  
  # Calculate survival probability
  sl_fit$density_failure_1$hazard_to_survival()
  sl_fit$density_failure_0$hazard_to_survival()
  
  k_grid <- 1:max(mydata_df$Time)
  sl_fit$density_failure_1$t <- k_grid
  sl_fit$density_failure_0$t <- k_grid
  
  sl_density_failure_1_marginal <- sl_fit$density_failure_1$clone(deep = TRUE)
  sl_density_failure_0_marginal <- sl_fit$density_failure_0$clone(deep = TRUE)
  # (sl_fit$density_failure_1$survival)[1,] 
  TSHEE_mydata_trt_1_survival_prob <- (sl_fit$density_failure_1$survival)
  TSHEE_mydata_trt_0_survival_prob <- (sl_fit$density_failure_0$survival)
  
  
  # Get the predicted median surviavl time from TSHEE model
  predict_TSHEE_mydata_median_survival_trt_1 <- rep(NA, dim(mydata$X)[1])
  predict_TSHEE_mydata_median_survival_trt_0 <- rep(NA, dim(mydata$X)[1])
  for (j in 1:dim(mydata$X)[1]){
    predict_TSHEE_mydata_median_survival_trt_1_index <- which(abs(TSHEE_mydata_trt_1_survival_prob[j,]-0.5) == min(abs(TSHEE_mydata_trt_1_survival_prob[j,]-0.5)))
    
    predict_TSHEE_mydata_median_survival_trt_1[j] <- k_grid[predict_TSHEE_mydata_median_survival_trt_1_index]
    
    predict_TSHEE_mydata_median_survival_trt_0_index <- which(abs(TSHEE_mydata_trt_0_survival_prob[j,]-0.5) == min(abs(TSHEE_mydata_trt_0_survival_prob[j,]-0.5)))
    
    predict_TSHEE_mydata_median_survival_trt_0[j] <- k_grid[predict_TSHEE_mydata_median_survival_trt_0_index]
  }
  
  # Test
  bias_ate <- tibble(p1_category = 1:50,
                     true_ate = rnorm(50, 14,0.5),
                     bias_ate = rnorm(50,0.25,0.05))
  
  # Calculate the bias
  true_ate <- tibble(true_mst_trt_2 = true_mst_trt_2, 
                     true_mst_trt_1 = true_mst_trt_1,
                     p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(true_ate = true_mst_trt_2 -  true_mst_trt_1) %>% 
    summarise(true_ate = mean(true_ate)) %>% 
    ungroup()
  
  TSHEE_ate <- tibble(predict_mst_trt_2 = predict_TSHEE_mydata_median_survival_trt_1, 
                        predict_mst_trt_1 = predict_TSHEE_mydata_median_survival_trt_0,
                        p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(ate = predict_mst_trt_2 -  predict_mst_trt_1) %>%
    # pull(ate) %>% mean
    summarise(ate = mean(ate)) 
  
  Bias_ate <- true_ate %>% 
    inner_join(TSHEE_ate) %>% 
    mutate(bias = true_ate - ate)
  return(bias_ate)
  
}