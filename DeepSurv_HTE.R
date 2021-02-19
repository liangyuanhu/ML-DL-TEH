for (j in 1:200){
  # For scenario 1: N = 5000, PH, 20% censor, heterogeneous setting (i), strong overlap-------------- 
  set.seed(j)
  mydata <- data_gen_censor(n=5000, p=10,  PH=TRUE, censor = "20%", setting = 1, overlap = "strong")
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
  # load the partial hazards result calculated in Python
  partial_trt_1 <- read_csv(paste0("partial_hazards_trt_1_simulation_",j,".csv"))
  partial_trt_0 <- read_csv(paste0("partial_hazards_trt_0_simulation_",j,".csv"))
  
  mydata_trt_1 <- mydata_trt_1 %>% 
    bind_cols(partial_trt_1)
  mydata_trt_0 <- mydata_trt_0 %>% 
    bind_cols(partial_trt_0)
  time_evaluated <- seq(1,50, 0.01)
  # Calcuate the survival probability for each subject at each time_evaluated
  deepsurv_mydata_trt_1_survival_prob <- matrix(nrow = dim(mydata_trt_1)[1], ncol = length(time_evaluated))
  deepsurv_mydata_trt_0_survival_prob <- matrix(nrow = dim(mydata_trt_0)[1], ncol = length(time_evaluated))
  for (m in length(time_evaluated)){
    basehaz_trt_1 <- gbm::basehaz.gbm(t = mydata_trt_1$Time, delta = mydata_trt_1$Event, f.x = log(-mydata_trt_1$partial_hazards), t.eval = time_evaluated[m])
    basehaz_trt_0 <- gbm::basehaz.gbm(t = mydata_trt_0$Time, delta = mydata_trt_0$Event, f.x = log(-mydata_trt_0$partial_hazards), t.eval = time_evaluated[m])
    deepsurv_mydata_trt_1_survival_prob[,m] <- exp(-(basehaz_trt_1 * exp(mydata_trt_1$partial_hazards)))
    deepsurv_mydata_trt_0_survival_prob[,m] <- exp(-(basehaz_trt_0 * exp(mydata_trt_0$partial_hazards)))
  }
  predict_deepsurv_mydata_median_survival_trt_1 <- rep(NA, dim(mydata_trt_1)[1])
  predict_deepsurv_mydata_median_survival_trt_0 <- rep(NA, dim(mydata_trt_0)[1])
  
  # Get the median survival time for each subject
  for (l in 1:dim(mydata_trt_1)[1]){
    predict_deepsurv_mydata_median_survival_trt_1_index <- which(abs(deepsurv_mydata_trt_1_survival_prob[l,]-0.5) == min(abs(deepsurv_mydata_trt_1_survival_prob[l,]-0.5)))
    predict_deepsurv_mydata_median_survival_trt_1[j] <- time_evaluated[predict_deepsurv_mydata_median_survival_trt_1_index]
  }
  for (l in 1:dim(mydata_trt_0)[1]){
    predict_deepsurv_mydata_median_survival_trt_0_index <- which(abs(deepsurv_mydata_trt_0_survival_prob[l,]-0.5) == min(abs(deepsurv_mydata_trt_0_survival_prob[l,]-0.5)))
    predict_deepsurv_mydata_median_survival_trt_0[j] <- time_evaluated[predict_deepsurv_mydata_median_survival_trt_0_index]
  }
  # Test
  bias_ate <- tibble(p1_category = 1:50,
                     true_ate = rnorm(50, 12.3,0.5),
                     bias_ate = rnorm(50,0.21,0.05),
                     regret = case_when(as.numeric(rbernoulli(50, 0.5)) == 1 ~ rnorm(50, 0.25, 0.04),
                                        TRUE ~ 0))
  
  # Calculate the bias and regret
  true_ate <- tibble(true_mst_trt_2 = true_mst_trt_2, 
                     true_mst_trt_1 = true_mst_trt_1,
                     p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(true_ate = true_mst_trt_2 -  true_mst_trt_1) %>% 
    summarise(true_ate = mean(true_ate)) %>% 
    ungroup()
  
  deepsurv_ate <- tibble(predict_mst_trt_2 = predict_deepsurv_mydata_median_survival_trt_1, 
                        predict_mst_trt_1 = predict_deepsurv_mydata_median_survival_trt_0,
                        p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(ate = predict_mst_trt_2 -  predict_mst_trt_1) %>%
    # pull(ate) %>% mean
    summarise(ate = mean(ate)) 
  
  Bias_ate <- true_ate %>% 
    inner_join(deepsurv_ate) %>% 
    mutate(bias = true_ate - ate,
           regret = ifelse(sign(true_ate,) == sign(ate), 0, bias))
  return(bias_ate)
  
}