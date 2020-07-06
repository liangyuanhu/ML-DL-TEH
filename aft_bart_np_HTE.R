aft_bart_np_HTE <- function(){
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
  
  # Fit AFTrees model 
  AFTrees_mydata_trt_1 = AFTrees(
    x.train = mydata_trt_1 %>% select(-Time, -Event) %>% as.matrix(),
    y.train =  mydata_trt_1$Time,
    status = mydata_trt_1$Event,
    nskip = 100,
    ndpost = 1000,
    x.test = mydata_trt_0 %>% select(-Time, -Event) %>% as.matrix(),
    nonparametric = T
  )
  
  AFTrees_mydata_trt_0 = AFTrees(
    x.train = mydata_trt_0 %>% select(-Time, -Event) %>% as.matrix(),
    y.train =  mydata_trt_0$Time,
    status = mydata_trt_0$Event,
    nskip = 100,
    ndpost = 1000,
    x.test = mydata_trt_1 %>% select(-Time, -Event) %>% as.matrix()
  )
  
  # Calculate survival probability
  AFTrees_mydata_trt_1_trt_1_survival_prob <- AFTrees_SurvivalProb(object = AFTrees_mydata_trt_1, train.only = T)
  AFTrees_mydata_trt_1_trt_0_survival_prob <- AFTrees_SurvivalProb(object = AFTrees_mydata_trt_1, test.only =T) 
  AFTrees_mydata_trt_1_survival_prob <- rbind(AFTrees_mydata_trt_1_trt_1_survival_prob$Surv.train, AFTrees_mydata_trt_1_trt_0_survival_prob$Surv.test)
  AFTrees_mydata_trt_0_trt_0_survival_prob <- AFTrees_SurvivalProb(object = AFTrees_mydata_trt_0, train.only = T)
  AFTrees_mydata_trt_0_trt_1_survival_prob <- AFTrees_SurvivalProb(object = AFTrees_mydata_trt_0, test.only = T)
  AFTrees_mydata_trt_0_survival_prob <- rbind(AFTrees_mydata_trt_0_trt_1_survival_prob$Surv.test, AFTrees_mydata_trt_0_trt_0_survival_prob$Surv.train)
  
  mydata_trt_1_time <- AFTrees_mydata_trt_1_trt_1_survival_prob$time.points
  mydata_trt_0_time <- AFTrees_mydata_trt_0_trt_0_survival_prob$time.points
  # Get the predicted median surviavl time from AFTrees model
  predict_AFTrees_mydata_median_survival_trt_1 <- rep(NA, dim(mydata$X)[1])
  predict_AFTrees_mydata_median_survival_trt_0 <- rep(NA, dim(mydata$X)[1])
  for (j in 1:dim(mydata$X)[1]){
    predict_AFTrees_mydata_median_survival_trt_1_index <- which(abs(AFTrees_mydata_trt_1_survival_prob[j,]-0.5) == min(abs(AFTrees_mydata_trt_1_survival_prob[j,]-0.5)))
    
    predict_AFTrees_mydata_median_survival_trt_1[j] <- mydata_trt_1_time[predict_AFTrees_mydata_median_survival_trt_1_index]
    
    predict_AFTrees_mydata_median_survival_trt_0_index <- which(abs(AFTrees_mydata_trt_0_survival_prob[j,]-0.5) == min(abs(AFTrees_mydata_trt_0_survival_prob[j,]-0.5)))
    
    predict_AFTrees_mydata_median_survival_trt_0[j] <- mydata_trt_0_time[predict_AFTrees_mydata_median_survival_trt_0_index]
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
  
  AFTrees_ate <- tibble(predict_mst_trt_2 = predict_AFTrees_mydata_median_survival_trt_1, 
                        predict_mst_trt_1 = predict_AFTrees_mydata_median_survival_trt_0,
                        p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(ate = predict_mst_trt_2 -  predict_mst_trt_1) %>%
    # pull(ate) %>% mean
    summarise(ate = mean(ate)) 
  
  Bias_ate <- true_ate %>% 
    inner_join(AFTrees_ate) %>% 
    mutate(bias = true_ate - ate)
  return(bias_ate)
  
}
  
  
  

