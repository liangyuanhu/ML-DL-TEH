aft_lognormal_HTE <- function(){
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
  # Data for treatment group
  mydata_trt_1 <- mydata$X %>% as_tibble() %>%  
    filter(W == 1) %>% 
    select(-W) %>% 
    as.data.frame()
  # Data for control group
  mydata_trt_0 <- mydata$X %>% as_tibble() %>%  
    filter(W == 0) %>% 
    select(-W) %>% 
    as.data.frame()
  
  mydata_aft_log_normal_linear_trt_1 <- survreg(Surv(Time, Event) ~ ., dist='lognormal', data = mydata_trt_1)
  mydata_aft_log_normal_linear_trt_0 <- survreg(Surv(Time, Event) ~ ., dist='lognormal', data = mydata_trt_0)
  
  # Predict the survival probability if everyone received treatemnt = 1
  predict_mydata_aft_log_normal_linear_trt_1_trt_1 <- predict(mydata_aft_log_normal_linear_trt_1, newdata = mydata_trt_1, type='quantile', p =  (1:98/100), se=TRUE) 
  predict_mydata_aft_log_normal_linear_trt_1_trt_0 <- predict(mydata_aft_log_normal_linear_trt_1, newdata = mydata_trt_0, type='quantile', p =  (1:98/100), se=TRUE) 
  predict_mydata_aft_log_normal_linear_trt_1 <- predict_mydata_aft_log_normal_linear_trt_1_trt_1$fit %>% 
    rbind(predict_mydata_aft_log_normal_linear_trt_1_trt_0$fit) 
  
  # Predict the survival probability if everyone received treatemnt = 0
  predict_mydata_aft_log_normal_linear_trt_0_trt_1 <- predict(mydata_aft_log_normal_linear_trt_0, newdata = mydata_trt_1, type='quantile', p =  (1:98/100), se=TRUE) 
  predict_mydata_aft_log_normal_linear_trt_0_trt_0 <- predict(mydata_aft_log_normal_linear_trt_0, newdata = mydata_trt_0, type='quantile', p =  (1:98/100), se=TRUE) 
  predict_mydata_aft_log_normal_linear_trt_0 <- predict_mydata_aft_log_normal_linear_trt_0_trt_1$fit %>% 
    rbind(predict_mydata_aft_log_normal_linear_trt_0_trt_0$fit) 
  
  # Get the predicted median surviavl time from AFT model
  predict_aft_log_normal_linear_mydata_median_survival_trt_1 <- rep(NA, dim(mydata$X)[1])
  predict_aft_log_normal_linear_mydata_median_survival_trt_0 <- rep(NA, dim(mydata$X)[1])
  for (j in 1:dim(mydata$X)[1]){
    predict_aft_log_normal_linear_mydata_median_survival_trt_1[j] <- predict_mydata_aft_log_normal_linear_trt_1[j,][which((1:98/100) == 0.5)]
    
    predict_aft_log_normal_linear_mydata_median_survival_trt_0[j] <- predict_mydata_aft_log_normal_linear_trt_0[j,][which((1:98/100) == 0.5)]
    
  }
  
  # Test
  bias_ate <- tibble(p1_category = 1:50,
         true_ate = rnorm(50, 14.5,0.5),
         bias_ate = rnorm(50,3.25,0.5))
  
  # Calculate the bias
  true_ate <- tibble(true_mst_trt_2 = true_mst_trt_2, 
                     true_mst_trt_1 = true_mst_trt_1,
                     p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(true_ate = true_mst_trt_2 -  true_mst_trt_1) %>% 
    summarise(true_ate = mean(true_ate)) %>% 
    ungroup()
  
  aft_log_normal_linear_ate <- tibble(predict_mst_trt_2 = predict_aft_log_normal_linear_mydata_median_survival_trt_1, 
                                      predict_mst_trt_1 = predict_aft_log_normal_linear_mydata_median_survival_trt_0,
                                      p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(ate = predict_mst_trt_2 -  predict_mst_trt_1) %>%
    # pull(ate) %>% mean
    summarise(ate = mean(ate)) 
  
  Bias_ate <- true_ate %>% 
    inner_join(aft_log_normal_linear_ate) %>% 
    mutate(bias = true_ate - ate)
  
  return(bias_ate)
}




