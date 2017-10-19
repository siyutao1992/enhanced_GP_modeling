getwd() # Please run the script at the directory of this file
# to do that, exit Rstudio and then click on this file in the relevant folder, then...
# the working directory is automatically set in this folder

library(R.matlab)
library(pracma)
library(lhs)
library(randtoolbox)

library(GPfit)
library(GPM)

setwd("./GPM_1.0_source/")
source("CovMat.R")
source("NLogL.R")
setwd("../")

# Start running for examples 1-7-----------------------------------------------

setwd("./datasets&results/")

NAMEs = c('Ex_1_10d', 'Ex_1_20d', 'Ex_1_40d',
          'Ex_2_10d', 'Ex_2_20d', 'Ex_2_40d',
          'Ex_3_10d', 'Ex_3_20d', 'Ex_3_40d',
          'Ex_4_10d', 'Ex_4_20d', 'Ex_4_40d')#,
          # 'Ex_5_10d', 'Ex_5_20d', 'Ex_5_40d',
          # 'Ex_6_10d', 'Ex_6_20d', 'Ex_6_40d',
          # 'Ex_7_10d', 'Ex_7_20d', 'Ex_7_40d') # array of names of datasets (not include '.mat')

N_replicates = 20 # please set the desired num of reps

GPfit_run = FALSE # switch for running GPfit

GPM_run = TRUE # switch for running GPM

# ==================== parameters all set now =========================

rep_ind_set = 1:N_replicates

len_NAMEs = length(NAMEs)

for (name_ind in 1:len_NAMEs) {

  NAME = NAMEs[name_ind]

  temp = readMat(strcat(NAME, '.mat'));

  X = get('temp')[[1]][[1]]
  Y = get('temp')[[1]][[2]]
  X_test = as.matrix(get('temp')[[1]][[3]])
  Y_test = as.matrix(get('temp')[[1]][[4]])
  Min = as.matrix(get('temp')[[1]][[5]])
  Max = as.matrix(get('temp')[[1]][[6]])

  GPfit_Err = matrix(0, nrow = N_replicates, ncol = 7)
  GPfit_Cost = matrix(0, nrow = N_replicates, ncol = 1)
  GPfit_All_models = list()

  GPfit_nugget = matrix(0, nrow = N_replicates, ncol = 1)
  GPfit_obj_func_value_conv = matrix(0, nrow = N_replicates, ncol = 1)

  GPM_Err = matrix(0, nrow = N_replicates, ncol = 7)
  GPM_Cost = matrix(0, nrow = N_replicates, ncol = 1)
  GPM_All_models = list()

  GPM_nugget = matrix(0, nrow = N_replicates, ncol = 1)
  GPM_obj_func_value = matrix(0, nrow = N_replicates, ncol = 1)

  for (i in rep_ind_set){
    X_i = as.matrix(X[, , i])
    X_iN = scale(X_i, Min, Max - Min);
    Y_i = as.matrix(Y[, , i])

    # GPfit part start --------------------------------------------------------------
    if (GPfit_run) {
      # fit the GPfit model
      ptm = proc.time()
      GPfit_Model = GP_fit(X_iN, Y_i, corr = list(type="exponential",power=2))
      temp = proc.time() - ptm
      GPfit_Cost[i] = temp["elapsed"]
      cat("\n","GPfit: rep#", i, "-fitted with time cost",GPfit_Cost[i], "(sec)", "\n" )
  
      # evaluate the prediction accuracy
      X_testN = scale(X_test, Min, Max - Min);
      temp = predict.GP(GPfit_Model,X_testN)
  
      # mean squared error
      GPfit_Err[i, 1] = mean((temp$Y_hat - Y_test)**2)
      # std of squared error
      GPfit_Err[i, 2] = std((temp$Y_hat - Y_test)**2)
      # maximum absolute error
      GPfit_Err[i, 3] = max(abs(temp$Y_hat - Y_test))
      # mean absolute error
      GPfit_Err[i, 4] = mean(abs(temp$Y_hat - Y_test))
      # std absolute error
      GPfit_Err[i, 5] = std(abs(temp$Y_hat - Y_test))
      # negatively oriented mean interval score for alpha = 0.05
      # see: Gneiting and Raftery, 2007
      part1 = temp$Y_hat - Y_test - 1.96*sqrt(temp$MSE)
      part2 = Y_test - temp$Y_hat - 1.96*sqrt(temp$MSE)
      GPfit_Err[i, 6] = mean(3.92*sqrt(temp$MSE) + 40*part1*(part1>0) + 40*part2*(part2>0))
      GPfit_Err[i, 7] = std(3.92*sqrt(temp$MSE) + 40*part1*(part1>0) + 40*part2*(part2>0))
  
      # get nugget and obj func value
      w = as.vector(GPfit_Model$beta)
      k = nrow(X_iN)
      q = ncol(Y_i)
      Fn = matrix(data = 1, nrow = k, ncol = 1)
      Ymin = apply(Y_i, 2, min)
      Yrange = apply(Y_i, 2, max) - Ymin
      Y_iN = t((t(Y_i)-Ymin)/Yrange)
  
      CovType='G'
      R = CovMat(X_iN, X_iN, CovType, w)
      Raw_MinEig = sort(eigen(R, symmetric = TRUE, only.values = TRUE)$values)[1]
      nugget = GPfit_Model$delta
      MinEig = nugget+Raw_MinEig
  
      GPfit_nugget[i] = nugget
      GPfit_obj_func_value_conv[i] = NLogL(w, X_iN, Y_iN, CovType, MinEig, Fn, k, q)
  
      tmp = strcat("GPfit_Model_", num2str(i, 0))
      GPfit_All_models[[tmp]] = GPfit_Model
  
  
      save(GPfit_All_models, GPfit_Cost, GPfit_Err, GPfit_obj_func_value_conv, GPfit_nugget,
           file = paste("GPfit_all_data_", NAME, "_reps1thru", toString(i), ".RData", sep = ""))
      if(i>1){
        file.remove(paste("GPfit_all_data_", NAME, "_reps1thru", toString(i-1), ".RData", sep = ""))
      }
    }
    
    # GPfit part end --------------------------------------------------------------

    # GPM part start --------------------------------------------------------------

    if(GPM_run) {
      # fit the GPM model
      GPM_Model = Fit(X_iN, Y_i) # "Progress" controls whether print out
      GPM_Cost[i] = GPM_Model$Details$Cost
      cat("\n","GPM: rep#", i, "-fitted with time cost",GPM_Cost[i], "(sec)","\n" )
  
      # evaluate the prediction accuracy
      X_testN = scale(X_test, Min, Max - Min);
      temp = Predict(X_testN, GPM_Model, MSE_on = 1)
  
      # mean squared error
      GPM_Err[i, 1] = mean((temp$YF - Y_test)**2)
      # std of squared error
      GPM_Err[i, 2] = std((temp$YF - Y_test)**2)
      # maximum absolute error
      GPM_Err[i, 3] = max(abs(temp$YF - Y_test))
      # mean absolute error
      GPM_Err[i, 4] = mean(abs(temp$YF - Y_test))
      # std absolute error
      GPM_Err[i, 5] = std(abs(temp$YF - Y_test))
      # negatively oriented mean interval score for alpha = 0.05
      # see: Gneiting and Raftery, 2007
      part1 = temp$YF - Y_test - 1.96*sqrt(temp$MSE)
      part2 = Y_test - temp$YF - 1.96*sqrt(temp$MSE)
      GPM_Err[i, 6] = mean(3.92*sqrt(temp$MSE) + 40*part1*(part1>0) + 40*part2*(part2>0))
      GPM_Err[i, 7] = std(3.92*sqrt(temp$MSE) + 40*part1*(part1>0) + 40*part2*(part2>0))
  
      GPM_nugget[i] = GPM_Model$Details$Nug_opt
      GPM_obj_func_value[i] = GPM_Model$Details$MinNLogL
  
      tmp = strcat("GPM_Model_", num2str(i, 0))
      GPM_All_models[[tmp]] = GPM_Model
  
      save(GPM_All_models, GPM_Cost, GPM_Err, GPM_obj_func_value, GPM_nugget,
           file = paste("GPM_all_data_", NAME, "_reps1thru", toString(i), ".RData", sep = ""))
      if(i>1){
        file.remove(paste("GPM_all_data_", NAME, "_reps1thru", toString(i-1), ".RData", sep = ""))
      }
    }
    
    # GPM part end --------------------------------------------------------------

  }
  # print out the results
  if(GPM_run) {
    cat('\n------------------------results-------------------------\n',
        NAME, '\n------------------------results-------------------------\n')
    cat(mean(GPM_Cost),'\n')
    cat(mean(GPM_Err[,1]),'\n')
    cat(mean(GPM_Err[,3]),'\n')
    cat(mean(GPM_Err[,6]),'\n')
  }
}

setwd("../") # back to the initial wd

# End running for examples 1-7-----------------------------------------------
