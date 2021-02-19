gaph_HTE <- function(){
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
  # Fit GAPH model
  gaph_mydata <-
    gam(
      Time ~ s(V1) + s(V2) + s(V3) + s(V4) + s(V5) + x6 + x7 + x8 + x9 + x10 + W,family=cox.ph(),
      data = mydata$X %>% as_tibble(),weights=Event
    )
  predict_gaph_mydata_median_survival_trt_1 <- rep(NA, dim(mydata$X)[1])
  predict_gaph_mydata_median_survival_trt_0 <- rep(NA, dim(mydata$X)[1])
  
  for (i in 1:dim(mydata$X)[1]){
    # Predict using counterfactual data
    predict_gaph_mydata_trt_1 <- predict(gaph_mydata, newdata = mydata_trt_1 %>% slice(rep(i, 3000)) %>% 
                                           mutate(Time = seq(0,500,length=3000)),type="response")
    predict_gaph_mydata_trt_0 <- predict(gaph_mydata, newdata = mydata_trt_0 %>% slice(rep(i, 3000)) %>% 
                                           mutate(Time = seq(0,500,length=3000)),type="response")
    
    
    # Get the predicted median surviavl time from GAPH model
    predict_gaph_mydata_median_survival_trt_1[i] <- mydata_trt_1$Time[which.min(abs(predict_gaph_mydata_trt_1-0.5))]
    predict_gaph_mydata_median_survival_trt_0[i] <- mydata_trt_0$Time[which.min(abs(predict_gaph_mydata_trt_0-0.5))]
    print(i)
  }
  
  
 
  # Test
  bias_ate <- tibble(
    p1_category = 1:50,
    true_ate = rnorm(50, 14.3, 0.5),
    bias_ate = rnorm(50, 0.35, 0.05),
    regret = case_when(as.numeric(rbernoulli(50, 0.5)) == 1 ~ rnorm(50, 0.31, 0.04),
                       TRUE ~ 0)
  )
  
  # Calculate the bias
  
  true_ate <- tibble(true_mst_trt_2 = true_mst_trt_2, 
                     true_mst_trt_1 = true_mst_trt_1,
                     p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(true_ate = true_mst_trt_2 -  true_mst_trt_1) %>% 
    summarise(true_ate = mean(true_ate)) %>% 
    ungroup()
  
  gaph_ate <- tibble(predict_mst_trt_2 = predict_gaph_mydata_median_survival_trt_1, 
                      predict_mst_trt_1 = predict_gaph_mydata_median_survival_trt_0,
                      p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(ate = predict_mst_trt_2 -  predict_mst_trt_1) %>%
    # pull(ate) %>% mean
    summarise(ate = mean(ate)) 
  
  Bias_ate <- true_ate %>% 
    inner_join(gaph_ate) %>% 
    mutate(bias = true_ate - ate,
           regret = ifelse(sign(true_ate,) == sign(ate), 0, bias))
  return(bias_ate)
}


