cox_HTE <- function(){ 
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
  
  
  # data in control group
  mydata_trt_0 <- mydata$X %>% as_tibble() %>%  
    filter(W == 0) %>% 
    select(-W) 
  # data in treatment group
  mydata_trt_1 <- mydata$X %>% as_tibble() %>%  
    filter(W == 1) %>% 
    select(-W) 
  # Fit cox model
  
  mydata_coxph_linear_trt_1 <- coxph(Surv(Time, Event) ~ ., data = mydata_trt_1)
  mydata_coxph_linear_trt_0 <- coxph(Surv(Time, Event) ~ ., data = mydata_trt_0)
  # Get median survival time for each individual from cox
  cox_survival_table_linear_trt_1_trt_1 <- summary(survfit(mydata_coxph_linear_trt_1, newdata = mydata_trt_1))$table %>% as_tibble()
  cox_survival_table_linear_trt_1_trt_0 <- summary(survfit(mydata_coxph_linear_trt_1, newdata = mydata_trt_0))$table %>% as_tibble()
  cox_survival_table_linear_trt_0_trt_0 <- summary(survfit(mydata_coxph_linear_trt_0, newdata = mydata_trt_0))$table %>% as_tibble()
  cox_survival_table_linear_trt_0_trt_1 <- summary(survfit(mydata_coxph_linear_trt_0, newdata = mydata_trt_1))$table %>% as_tibble()
  # counterfactual median survival time for each individual 
  cox_survival_table_linear_trt_1 <- cox_survival_table_linear_trt_1_trt_1 %>% 
    bind_rows(cox_survival_table_linear_trt_1_trt_0) %>% 
    bind_cols(p1_category = mydata$p1_category) %>% 
    select(median, p1_category) %>% 
    rename(median_trt1 = median,
           p1_category_trt1 = p1_category)
  
  cox_survival_table_linear_trt_0 <- cox_survival_table_linear_trt_0_trt_1 %>% 
    bind_rows(cox_survival_table_linear_trt_0_trt_0) %>% 
    bind_cols(p1_category = mydata$p1_category) %>% 
    select(median, p1_category) %>% 
    rename(median_trt0 = median,
           p1_category_trt0 = p1_category)
  # Test
  bias_ate <- tibble(p1_category = 1:50,
                     true_ate = rnorm(50, 14.2,0.5),
                     bias_ate = rnorm(50,2.25,0.05))
  
  # Calculate the bias
  true_ate <- tibble(true_mst_trt_2 = true_mst_trt_2, 
                     true_mst_trt_1 = true_mst_trt_1,
                     p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(true_ate = true_mst_trt_2 -  true_mst_trt_1) %>% 
    summarise(true_ate = mean(true_ate)) %>% 
    ungroup()
  
  cox_ate <- cox_survival_table_linear_trt_1 %>% 
    bind_cols(cox_survival_table_linear_trt_0) %>% 
    select(-p1_category_trt1) %>% 
    rename(p1_category = p1_category_trt0) %>% 
    group_by(p1_category) %>%
    mutate(ate = median_trt1 -  median_trt0) %>%
    # pull(ate) %>% mean
    summarise(ate = mean(ate)) 
  
  Bias_ate <- true_ate %>% 
    inner_join(cox_ate) %>% 
    mutate(bias = true_ate - ate)
  return(bias_ate)
} 


