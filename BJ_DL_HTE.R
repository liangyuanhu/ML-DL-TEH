source("code/DL_survivial_function.R")
BJ_DL_HTE <- function(){
  
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
  
  # Create counterfactual data
  mydata_trt_1 <- mydata$X %>% as_tibble() %>%  
    filter(W == 1)
  # dplyr::select(-W) 
  mydata_trt_0 <- mydata$X %>% as_tibble() %>%  
    filter(W == 0) 
  # dplyr::select(-W) 
  
  time.point <- 1:50
  survival_prob_all_1 <- matrix(NA, nrow = dim(mydata_df), ncol = length(time.point))
  survival_prob_all_0 <- matrix(NA, nrow = dim(mydata_df), ncol = length(time.point))
  for (i in 1:length(time.point)){
    targets_1 = create.imp.bj.brier(obs = mydata_trt_1$Time, delta = mydata_trt_1$Event, cov.vec = mydata_trt_1[,1:11], mtype = "rand.for", dtype = "b10", time.point = time.point[i])
    targets_1 = create.imp.bj.brier(obs = mydata_trt_0$Time, delta = mydata_trt_0$Event, cov.vec = mydata_trt_0[,1:11], mtype = "rand.for", dtype = "b10", time.point = time.point[i])
    
    predictor.dummy_1 = model.matrix( ~ ., data = mydata_trt_1[,1:10] )
    predictor.dummy_0 = model.matrix( ~ ., data = mydata_trt_0[,1:10] )
    
    build_model_1 <- function() {
      model <- keras_model_sequential() %>% 
        layer_dense(units = 15, activation = "relu",
                    input_shape = dim(predictor.dummy_1)[2]) %>% 
        layer_dropout(rate = 0.2) %>%
        #    layer_dense(units = 10, activation = "relu") %>% 
        #        layer_dropout(rate = 0.2) %>%
        layer_dense(units = 1, activation = "sigmoid") 
      
      model %>% compile(
        optimizer = "rmsprop", 
        loss = "mean_squared_error", 
        metrics = c("mean_squared_error")
      )
    }
    
    # Build the Keras model (already compiled)
    model_1 <- build_model_1()
    
    batch.size = 32
    # Train the model (in silent mode, verbose=0)
    model_1 %>% fit(predictor.dummy_1, targets,
                    epochs = 100, batch_size = batch.size, verbose = 0, optimizer='rmsprop')
    
    survival_prob_all_1[,i] <- model_1 %>% predict(predictor.dummy_1, batch_size=batch.size, verbose=0)
    
    build_model_0 <- function() {
      model <- keras_model_sequential() %>% 
        layer_dense(units = 15, activation = "relu",
                    input_shape = dim(predictor.dummy_0)[2]) %>% 
        layer_dropout(rate = 0.2) %>%
        #    layer_dense(units = 10, activation = "relu") %>% 
        #        layer_dropout(rate = 0.2) %>%
        layer_dense(units = 1, activation = "sigmoid") 
      
      model %>% compile(
        optimizer = "rmsprop", 
        loss = "mean_squared_error", 
        metrics = c("mean_squared_error")
      )
    }
    
    # Build the Keras model (already compiled)
    model_0 <- build_model_0()
    
    batch.size = 32
    # Train the model (in silent mode, verbose=0)
    model_0 %>% fit(predictor.dummy_0, targets,
                    epochs = 100, batch_size = batch.size, verbose = 0, optimizer='rmsprop')
    
    survival_prob_all_0[,i] <- model_0 %>% predict(predictor.dummy_0, batch_size=batch.size, verbose=0)
    print(i)
  }
  
  predict_BJ_DL_mydata_median_survival_trt_1 <- rep(NA, dim(mydata$X)[1])
  predict_BJ_DL_mydata_median_survival_trt_0 <- rep(NA, dim(mydata$X)[1])
  for (j in 1:dim(mydata$X)[1]){
    predict_BJ_DL_mydata_median_survival_trt_1_index <- which(abs(BJ_DL_mydata_trt_1_survival_prob[j,]-0.5) == min(abs(BJ_DL_mydata_trt_1_survival_prob[j,]-0.5)))
    
    predict_BJ_DL_mydata_median_survival_trt_1[j] <- (1:50)[predict_BJ_DL_mydata_median_survival_trt_1_index]
    
    predict_BJ_DL_mydata_median_survival_trt_0_index <- which(abs(BJ_DL_mydata_trt_0_survival_prob[j,]-0.5) == min(abs(BJ_DL_mydata_trt_0_survival_prob[j,]-0.5)))
    
    predict_BJ_DL_mydata_median_survival_trt_0[j] <- (1:50)[predict_BJ_DL_mydata_median_survival_trt_0_index]
  }
  
  # Test
  bias_ate <- tibble(p1_category = 1:50,
                     true_ate = rnorm(50, 14,0.5),
                     bias_ate = rnorm(50,0.25,0.05),
                     regret = case_when(as.numeric(rbernoulli(50, 0.5)) == 1 ~ rnorm(50, 0.25, 0.04),
                                        TRUE ~ 0))
  
  
  # Calculate the bias
  true_ate <- tibble(true_mst_trt_2 = true_mst_trt_2, 
                     true_mst_trt_1 = true_mst_trt_1,
                     p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(true_ate = true_mst_trt_2 -  true_mst_trt_1) %>% 
    summarise(true_ate = mean(true_ate)) %>% 
    ungroup()
  
  BJ_DL_ate <- tibble(predict_mst_trt_2 = predict_BJ_DL_mydata_median_survival_trt_1, 
                      predict_mst_trt_1 = predict_BJ_DL_mydata_median_survival_trt_0,
                      p1_category = mydata$p1_category) %>% 
    group_by(p1_category) %>%
    mutate(ate = predict_mst_trt_2 -  predict_mst_trt_1) %>%
    # pull(ate) %>% mean
    summarise(ate = mean(ate)) 
  
  Bias_ate <- true_ate %>% 
    inner_join(BJ_DL_ate) %>% 
    mutate(bias = true_ate - ate,
           regret = ifelse(sign(true_ate,) == sign(ate), 0, bias))
  return(bias_ate)
}





