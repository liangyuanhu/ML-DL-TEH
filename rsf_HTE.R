rsf_HTE <- function(){
  # True median survival time
  true_mst_trt_1 <- rep(NA, dim(mydata$X)[1])
  for (i in 1:dim(mydata$X)[1]){
    fun_trt_1 <-
      function(x) {
        exp(-(1000/mydata$lambda1 * exp(mydata$LP1[i]) * x) ^ mydata$eta_1)-0.5
      }
    true_mst_trt_1[i] <- uniroot(fun_trt_1, interval = c(0, dim(mydata$X)[1]))$root
    # print(i)
  }
  
  true_mst_trt_2 <- rep(NA, dim(mydata$X)[1])
  for (i in 1:dim(mydata$X)[1]){
    fun_trt_2 <-
      function(x) {
        exp(-(1000/mydata$lambda2 * exp(mydata$LP2[i]) * x) ^ mydata$eta_2)-0.5
      }
    true_mst_trt_2[i] <- uniroot(fun_trt_2, interval = c(0, dim(mydata$X)[1]))$root
    # print(i)
  }
  
  # Create counterfactual data
  mydata_trt_1 <- mydata$X %>% as_tibble() %>%  
    mutate(W = 1)
  
  mydata_trt_0 <- mydata$X %>% as_tibble() %>%  
    mutate(W = 0)
  # Fit RFSRC model
  rfsrc_mydata <-
    rfsrc(
      Surv(Time, Event) ~ .,
      data = mydata$X %>% as_tibble()
    )
  # Predict using counterfactual data
  predict_rfsrc_mydata_trt_1 <- predict(rfsrc_mydata, newdata = mydata_trt_1)
  predict_rfsrc_mydata_trt_0 <- predict(rfsrc_mydata, newdata = mydata_trt_0)
  
  
  predict_rfsrc_mydata_median_survival_trt_1 <- rep(NA,dim(mydata$X)[1])
  predict_rfsrc_mydata_median_survival_trt_0 <- rep(NA,dim(mydata$X)[1])
  # Get the predicted median surviavl time from RFSRC model
  for (j in 1:dim(mydata$X)[1]){
    
    predict_rfsrc_mydata_median_survival_trt_1_index <- which(abs(predict_rfsrc_mydata_trt_1$survival[j,]-0.5) == min(abs(predict_rfsrc_mydata_trt_1$survival[j,]-0.5)))
    
    predict_rfsrc_mydata_median_survival_trt_1[j] <- predict_rfsrc_mydata_trt_1$time.interest[predict_rfsrc_mydata_median_survival_trt_1_index]
    predict_rfsrc_mydata_median_survival_trt_0_index <- which(abs(predict_rfsrc_mydata_trt_0$survival[j,]-0.5) == min(abs(predict_rfsrc_mydata_trt_0$survival[j,]-0.5)))
    
    predict_rfsrc_mydata_median_survival_trt_0[j] <- predict_rfsrc_mydata_trt_0$time.interest[predict_rfsrc_mydata_median_survival_trt_0_index]
    
  }
  
  # Test
  bias_ate <- tibble(p1_category = 1:50,
                     true_ate = rnorm(50, 14.3,0.5),
                     bias_ate = rnorm(50,0.35,0.05))
  
  # Calculate the bias
  true_ate <- tibble(true_mst_trt_2 = true_mst_trt_2, 
                     true_mst_trt_1 = true_mst_trt_1,
                     p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(true_ate = true_mst_trt_2 -  true_mst_trt_1) %>% 
    summarise(true_ate = mean(true_ate)) %>% 
    ungroup()
  
  rfsrc_ate <- tibble(predict_mst_trt_2 = predict_rfsrc_mydata_median_survival_trt_1, 
                      predict_mst_trt_1 = predict_rfsrc_mydata_median_survival_trt_0,
                      p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(ate = predict_mst_trt_2 -  predict_mst_trt_1) %>%
    # pull(ate) %>% mean
    summarise(ate = mean(ate)) 
  
  Bias_ate <- true_ate %>% 
    inner_join(rfsrc_ate) %>% 
    mutate(bias = true_ate - ate)
  return(bias_ate)
}



