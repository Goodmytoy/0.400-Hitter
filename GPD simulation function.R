GPD_sim_test <- function(data, threshold, top_n=NA, intercept = T, result_loc, n_bootstrap =1000, iter_unit = 100){
  # browser()
  
  if(is.numeric(top_n)==T){
    top_n_threshold <- aggregate(BA~yearID, data, function(x) sort(x,decreasing =T)[top_n+1])
    threshold_of_year<-rep(NA, nrow(data))
    for(i in 1:nrow(data)){
      threshold_of_year[i] <- top_n_threshold[which(data[i,"yearID"] == top_n_threshold$yearID),"BA"]
    }
    data <- cbind(data, threshold_of_year)
    data <- data[data$BA > data$threshold_of_year,]
    threshold <- data$threshold_of_year
  }else{
    data <- data[data$BA > threshold,]
    threshold <- rep(threshold, nrow(data))
  }
  
  
  #check the location for saving results.
  if(dir.exists(result_loc)==F) dir.create(result_loc)
  
  ##Define Functions
  #GPD mean
  GPD_mean <- function(mu, sigma, xi){
    return(mu + sigma/(1-xi))
  }
  
  
  ##Find the best model
  cat("FInding the best model for 3 cases(scale_shape, scale, shape)\n")
  #model list
  iSpline_list_scale_shape <- vector("list", 16)
  iSpline_list_scale <- vector("list", 16)
  iSpline_list_shape <- vector("list", 16)
  
  #mle matrix
  mle_mat_scale_shape <- matrix(NA, ncol = 4, nrow = 4, 
                                dimnames=list(c("1 degree","2 degree", "3 degree", "4 degree"), c("1 knot","2 knots", "3 knots", "4 knots")))
  mle_mat_scale <- matrix(NA, ncol = 4, nrow = 4, 
                          dimnames=list(c("1 degree","2 degree", "3 degree", "4 degree"), c("1 knot","2 knots", "3 knots", "4 knots")))
  mle_mat_shape <- matrix(NA, ncol = 4, nrow = 4, 
                          dimnames=list(c("1 degree","2 degree", "3 degree", "4 degree"), c("1 knot","2 knots", "3 knots", "4 knots")))  
  
  knots_list <- list(round(quantile(data$year,0.5)), 
                     c(round(quantile(data$year,0.33)), round(quantile(data$year,0.66))),
                     c(round(quantile(data$year,0.25)), round(quantile(data$year,0.5)), round(quantile(data$year,0.75))),
                     c(round(quantile(data$year,0.2)), round(quantile(data$year,0.4)), round(quantile(data$year,0.6)), round(quantile(data$year,0.8)))
  )
  
  
  count <- 0
  for(i in 1:4){
    for(j in 1:4){
      cat("  Check model with N_knots = ", i, ", degree = ", j, "\n")
      count <- count + 1
      n_knots <- i
      dgr <- j
      
      #fitting
      gpd_ispline_scale_shape <- fevd(data$BA, threshold = threshold, type = "GP", method = "MLE",
                                      shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[i]],degree=dgr,intercept=intercept)[-1,]-1,
                                      scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[i]],degree=dgr,intercept=intercept)[-1,]-1
      )
      
      gpd_ispline_scale <- fevd(data$BA, threshold = threshold, type = "GP", method = "MLE",
                                scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[i]],degree=dgr,intercept=intercept)[-1,]-1
      )
      
      gpd_ispline_shape <- fevd(data$BA, threshold = threshold, type = "GP", method = "MLE",
                                shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[i]],degree=dgr,intercept=intercept)[-1,]-1
      )
      
      iSpline_list_scale_shape[[count]] <- gpd_ispline_scale_shape
      iSpline_list_scale[[count]] <- gpd_ispline_scale
      iSpline_list_shape[[count]] <- gpd_ispline_shape
      
      names(iSpline_list_scale_shape[count]) <- paste0("knots = ", i, ", degree = ", j, "\n")
      names(iSpline_list_scale[count]) <- paste0("knots = ", i, ", degree = ", j, "\n")
      names(iSpline_list_shape[count]) <- paste0("knots = ", i, ", degree = ", j, "\n")
      
      #column : knots, row : degree
      mle_mat_scale_shape[j,i] <- gpd_ispline_scale_shape$results$value
      mle_mat_scale[j,i] <- gpd_ispline_scale$results$value
      mle_mat_shape[j,i] <- gpd_ispline_shape$results$value
      
      #obtain parameter
      pars_scale_shape <- findpars(iSpline_list_scale_shape[[count]])
      pars_scale <- findpars(iSpline_list_scale[[count]])
      pars_shape <- findpars(iSpline_list_shape[[count]])
      
      
      GPD_mean_vec_scale_shape <- GPD_mean(mu = threshold, sigma = pars_scale_shape$scale, xi = pars_scale_shape$shape)
      GPD_mean_vec_scale <- GPD_mean(mu = threshold, sigma = pars_scale$scale, xi = pars_scale$shape)
      GPD_mean_vec_shape <- GPD_mean(mu = threshold, sigma = pars_shape$scale, xi = pars_shape$shape)
      
      #plot
      if(dir.exists(paste0(result_loc, "/scale_shape_par_plot")) == F) dir.create(paste0(result_loc, "/scale_shape_par_plot"))
      if(dir.exists(paste0(result_loc, "/scale_par_plot")) == F) dir.create(paste0(result_loc, "/scale_par_plot"))
      if(dir.exists(paste0(result_loc, "/shape_par_plot")) == F) dir.create(paste0(result_loc, "/shape_par_plot"))
      
      x_range <- data$year
      
      par_plot <- function(pars, x_range, knots, degree, type){
        # browser()
        if(type == "scale_shape"){
          plot_loc <- paste0(result_loc, "/scale_shape_par_plot")
        }else if(type == "scale"){
          plot_loc <- paste0(result_loc, "/scale_par_plot")
        }else if(type == "shape"){
          plot_loc <- paste0(result_loc, "/shape_par_plot")
        }
        
        GPD_mean_vec <- GPD_mean(mu=threshold, sigma = pars$scale, xi = pars$shape)
        
        png(paste0(plot_loc,"/",i,"knots, ", j, "degree.png"), 500, 500)
        par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
        plot(x_range, pars$scale, xlab = "year", ylab = "scale", main = "year - scale")
        plot(x_range, pars$shape, xlab = "year", ylab = "shape", main = "year - shape")
        plot(x_range, GPD_mean_vec, xlab = "year", ylab = "GPD mean", main = "year - GPD mean")
        dev.off()
        par(mfrow=c(1,1))
      }
      
      par_plot(pars_scale_shape, x_range, knots=i, degree=j, type="scale_shape")
      par_plot(pars_scale, x_range, knots=i, degree=j, type="scale")
      par_plot(pars_shape, x_range, knots=i, degree=j, type="shape")
    }
  }
  
  #saving MLE matrix
  cat("Saving MLE matrix. \n")
  write.csv(mle_mat_scale_shape, paste0(result_loc,"/mle_mat_scale_shape.csv"), row.names = F)
  write.csv(mle_mat_scale, paste0(result_loc,"/mle_mat_scale.csv"), row.names = F)
  write.csv(mle_mat_shape, paste0(result_loc,"/mle_mat_shape.csv"), row.names = F)
  
  #select best models  
  best_idx_scale_shape <- which.min(mle_mat_scale_shape)
  best_idx_scale <- which.min(mle_mat_scale)
  best_idx_shape <- which.min(mle_mat_shape)
  
  cat(best_idx_scale_shape, "\n")
  cat(best_idx_scale, "\n")
  cat(best_idx_scale, "\n")
  
  num_scale_shape <- which(mle_mat_scale_shape == min(mle_mat_scale_shape), arr.ind = T)
  num_scale <- which(mle_mat_scale == min(mle_mat_scale), arr.ind = T)
  num_shape <- which(mle_mat_shape == min(mle_mat_shape), arr.ind = T)
  
  
  gpd_constant <- fevd(data$BA, threshold = threshold, type = "GP", method = "MLE")
  
  best_model_scale_shape <- fevd(data$BA, threshold = threshold, type = "GP", method = "MLE",
                                 shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale_shape[2]]],degree=num_scale_shape[1],intercept=intercept)[-1,]-1,
                                 scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale_shape[2]]],degree=num_scale_shape[1],intercept=intercept)[-1,]-1
  )
  best_model_scale <- fevd(data$BA, threshold = threshold, type = "GP", method = "MLE",
                           scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale[2]]],degree=num_scale[1],intercept=intercept)[-1,]-1
  )
  best_model_shape <- fevd(data$BA, threshold = threshold, type = "GP", method = "MLE",
                           shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_shape[2]]],degree=num_shape[1],intercept=intercept)[-1,]-1
  )
  
  pars_constant <- findpars(gpd_constant)
  pars_scale_shape <- findpars(best_model_scale_shape)
  pars_scale <- findpars(best_model_scale)
  pars_shape <- findpars(best_model_shape)

  
  cat("The best models are found\n\n")
  cat("scale_shape model : ", num_scale_shape[1], " degree", ", ", num_scale_shape[2], " knots", "\n")
  cat("scale model : ", num_scale[1], " degree", ", ", num_scale[2], " knots", "\n")
  cat("shape model : ", num_shape[1], " degree", ", ", num_shape[2], " knots", "\n\n")
  
  #output.txt
  cat("scale_shape model : ", num_scale_shape[1], " degree", ", ", num_scale_shape[2], " knots", "\n",file=paste0(result_loc,"/best_models.txt"))
  cat("scale model : ", num_scale[1], " degree", ", ", num_scale[2], " knots", "\n" ,append = T, file=paste0(result_loc,"/best_models.txt"))
  cat("shape model : ", num_shape[1], " degree", ", ", num_shape[2], " knots", "\n\n", append = T, file=paste0(result_loc,"/best_models.txt"))
  cat("\n\n##############################################################################\n", append = T, file=paste0(result_loc,"/best_models.txt"))
  cat("Constant Model Output\n", append = T, file=paste0(result_loc,"/best_models.txt"))
  capture.output(gpd_constant, append = T, file=paste0(result_loc,"/best_models.txt"))
  cat("\n\n##############################################################################\n", append = T, file=paste0(result_loc,"/best_models.txt"))
  cat("Scale_Shape Model Output\n", append = T, file=paste0(result_loc,"/best_models.txt"))
  capture.output(best_model_scale_shape, append = T, file=paste0(result_loc,"/best_models.txt"))
  cat("\n\n##############################################################################\n", append = T, file=paste0(result_loc,"/best_models.txt"))
  cat("Scale Model Output\n", append = T, file=paste0(result_loc,"/best_models.txt"))
  capture.output(best_model_scale, append = T, file=paste0(result_loc,"/best_models.txt"))
  cat("\n\n##############################################################################\n", append = T, file=paste0(result_loc,"/best_models.txt"))
  cat("Shape Model Output\n", append = T, file=paste0(result_loc,"/best_models.txt"))
  capture.output(best_model_shape, append = T, file=paste0(result_loc,"/best_models.txt"))
  cat("\n\n##############################################################################\n", append = T, file=paste0(result_loc,"/best_models.txt"))
  #################################################################
  #Bootstraping
  cat("############################################################\n")
  cat("Start Bootstrap\n")
  
  
  #log likelihood vectors
  constant_loglik_from_constant <- rep(NA, n_bootstrap)
  constant_loglik_from_scale_shape <- rep(NA, n_bootstrap)
  constant_loglik_from_scale <- rep(NA, n_bootstrap)
  constant_loglik_from_shape <- rep(NA, n_bootstrap)
  
  scale_shape_loglik_from_constant <- rep(NA, n_bootstrap)
  scale_shape_loglik_from_scale_shape <- rep(NA, n_bootstrap)
  scale_shape_loglik_from_scale <- rep(NA, n_bootstrap)
  scale_shape_loglik_from_shape <- rep(NA, n_bootstrap)
  
  scale_loglik_from_constant <- rep(NA, n_bootstrap)
  scale_loglik_from_scale_shape <- rep(NA, n_bootstrap)
  scale_loglik_from_scale <- rep(NA, n_bootstrap)
  scale_loglik_from_shape <- rep(NA, n_bootstrap)
  
  shape_loglik_from_constant <- rep(NA, n_bootstrap)
  shape_loglik_from_scale_shape <- rep(NA, n_bootstrap)
  shape_loglik_from_scale <- rep(NA, n_bootstrap)
  shape_loglik_from_shape <- rep(NA, n_bootstrap)
  
  cat("Set Constant model as Null. \n")
  # Set Constant model as Null
  for(i in 1:n_bootstrap){
    if(i%%iter_unit==0) cat("Bootstrap iter", i, "\n")
    #data generating
    BA_constant <- rep(NA, nrow(data))
    BA_scale_shape <- rep(NA, nrow(data))
    BA_scale <- rep(NA, nrow(data))
    BA_shape <- rep(NA, nrow(data))
    
    for(j in 1:nrow(data)){
      BA_constant[j] <- revd(n=1, scale=pars_constant$scale[j], shape=pars_constant$shape[j], threshold=threshold[j], type="GP")
      BA_scale_shape[j] <- revd(n=1, scale=pars_scale_shape$scale[j], shape=pars_scale_shape$shape[j], threshold=threshold[j], type="GP")
      BA_scale[j] <- revd(n=1, scale=pars_scale$scale[j], shape=pars_scale$shape[j], threshold=threshold[j], type="GP")
      BA_shape[j] <- revd(n=1, scale=pars_shape$scale[j], shape=pars_shape$shape[j], threshold=threshold[j], type="GP")
    }
    bootstrap_data <- data.frame(year = data$year[data$BA>threshold],
                                 BA_constant = BA_constant,
                                 BA_scale_shape = BA_scale_shape,
                                 BA_scale = BA_scale,
                                 BA_shape = BA_shape)
    # bootstrap_data <- data.frame(year = data$year[data$BA>threshold],
    #                              BA_constant = rextRemes(x=gpd_constant, sum(data$BA>threshold)),
    #                              BA_scale_shape = rextRemes(x=best_model_scale_shape, 1),
    #                              BA_scale = rextRemes(x=best_model_scale, 1),
    #                              BA_shape = rextRemes(x=best_model_shape, 1))
    
    #constant models
    bootstrap_gpd_constant1 <- fevd(bootstrap_data$BA_constant, threshold = threshold, type = "GP", method = "MLE")
    bootstrap_gpd_constant2 <- fevd(bootstrap_data$BA_scale_shape, threshold = threshold, type = "GP", method = "MLE")
    bootstrap_gpd_constant3 <- fevd(bootstrap_data$BA_scale, threshold = threshold, type = "GP", method = "MLE")
    bootstrap_gpd_constant4 <- fevd(bootstrap_data$BA_shape, threshold = threshold, type = "GP", method = "MLE")
    
    constant_loglik_from_constant[i] <- -bootstrap_gpd_constant1$result$value
    constant_loglik_from_scale_shape[i] <- -bootstrap_gpd_constant2$result$value
    constant_loglik_from_scale[i] <- -bootstrap_gpd_constant3$result$value
    constant_loglik_from_shape[i] <- -bootstrap_gpd_constant4$result$value
    
    
    #scale, shape models
    bootstrap_gpd_scale_shape1 <- fevd(bootstrap_data$BA_constant, threshold = threshold, type = "GP", method = "MLE",
                                  shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale_shape[2]]],degree=num_scale_shape[1],intercept=intercept)[-1,]-1,
                                  scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale_shape[2]]],degree=num_scale_shape[1],intercept=intercept)[-1,]-1)
    bootstrap_gpd_scale_shape2 <- fevd(bootstrap_data$BA_scale_shape, threshold = threshold, type = "GP", method = "MLE",
                                       shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale_shape[2]]],degree=num_scale_shape[1],intercept=intercept)[-1,]-1,
                                       scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale_shape[2]]],degree=num_scale_shape[1],intercept=intercept)[-1,]-1)
    bootstrap_gpd_scale_shape3 <- fevd(bootstrap_data$BA_scale, threshold = threshold, type = "GP", method = "MLE",
                                       shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale_shape[2]]],degree=num_scale_shape[1],intercept=intercept)[-1,]-1,
                                       scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale_shape[2]]],degree=num_scale_shape[1],intercept=intercept)[-1,]-1)
    bootstrap_gpd_scale_shape4 <- fevd(bootstrap_data$BA_shape, threshold = threshold, type = "GP", method = "MLE",
                                       shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale_shape[2]]],degree=num_scale_shape[1],intercept=intercept)[-1,]-1,
                                       scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale_shape[2]]],degree=num_scale_shape[1],intercept=intercept)[-1,]-1)
    
    scale_shape_loglik_from_constant[i] <- -bootstrap_gpd_scale_shape1$result$value
    scale_shape_loglik_from_scale_shape[i] <- -bootstrap_gpd_scale_shape2$result$value
    scale_shape_loglik_from_scale[i] <- -bootstrap_gpd_scale_shape3$result$value
    scale_shape_loglik_from_shape[i] <- -bootstrap_gpd_scale_shape4$result$value
    
    
    #scale models
    bootstrap_gpd_scale1 <- fevd(bootstrap_data$BA_constant, threshold = threshold, type = "GP", method = "MLE",
                            scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale[2]]],degree=num_scale[1],intercept=T)[-1,]-1)
    bootstrap_gpd_scale2 <- fevd(bootstrap_data$BA_scale_shape, threshold = threshold, type = "GP", method = "MLE",
                                scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale[2]]],degree=num_scale[1],intercept=T)[-1,]-1)
    bootstrap_gpd_scale3 <- fevd(bootstrap_data$BA_scale, threshold = threshold, type = "GP", method = "MLE",
                                scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale[2]]],degree=num_scale[1],intercept=T)[-1,]-1)
    bootstrap_gpd_scale4 <- fevd(bootstrap_data$BA_shape, threshold = threshold, type = "GP", method = "MLE",
                                scale.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_scale[2]]],degree=num_scale[1],intercept=T)[-1,]-1)
    
    scale_loglik_from_constant[i] <- -bootstrap_gpd_scale1$result$value
    scale_loglik_from_scale_shape[i] <- -bootstrap_gpd_scale2$result$value
    scale_loglik_from_scale[i] <- -bootstrap_gpd_scale3$result$value
    scale_loglik_from_shape[i] <- -bootstrap_gpd_scale4$result$value
    
    #shape models
    bootstrap_gpd_shape1 <- fevd(bootstrap_data$BA_constant, threshold = threshold, type = "GP", method = "MLE",
                            shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_shape[2]]],degree=num_shape[1],intercept=T)[-1,]-1)
    bootstrap_gpd_shape2 <- fevd(bootstrap_data$BA_scale_shape, threshold = threshold, type = "GP", method = "MLE",
                                shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_shape[2]]],degree=num_shape[1],intercept=T)[-1,]-1)
    bootstrap_gpd_shape3 <- fevd(bootstrap_data$BA_scale, threshold = threshold, type = "GP", method = "MLE",
                                shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_shape[2]]],degree=num_shape[1],intercept=T)[-1,]-1)
    bootstrap_gpd_shape4 <- fevd(bootstrap_data$BA_shape, threshold = threshold, type = "GP", method = "MLE",
                                shape.fun = ~ iSpline(c(1899,data$year),knots=knots_list[[num_shape[2]]],degree=num_shape[1],intercept=T)[-1,]-1)
    
    shape_loglik_from_constant[i] <- -bootstrap_gpd_shape1$result$value
    shape_loglik_from_scale_shape[i] <- -bootstrap_gpd_shape2$result$value
    shape_loglik_from_scale[i] <- -bootstrap_gpd_shape3$result$value
    shape_loglik_from_shape[i] <- -bootstrap_gpd_shape4$result$value
  }
  
  ##NUll distribution
  cat("Calculating Null distribution. \n")
  #constant(null)
  constant_null_scale_shape_dist <- -2*(constant_loglik_from_constant - scale_shape_loglik_from_constant)
  constant_null_scale_dist <- -2*(constant_loglik_from_constant - scale_loglik_from_constant)
  constant_null_shape_dist <- -2*(constant_loglik_from_constant - shape_loglik_from_constant)
  
  #scale_shape(null)
  scale_shape_null_constant_dist <- -2*(scale_shape_loglik_from_scale_shape - constant_loglik_from_scale_shape)
  scale_shape_null_scale_dist <- -2*(scale_shape_loglik_from_scale_shape - scale_loglik_from_scale_shape)
  scale_shape_null_shape_dist <- -2*(scale_shape_loglik_from_scale_shape - shape_loglik_from_scale_shape)
  
  #scale(null)
  scale_null_constant_dist <- -2*(scale_loglik_from_scale - constant_loglik_from_scale)
  scale_null_scale_shape_dist <- -2*(scale_loglik_from_scale - scale_shape_loglik_from_scale)
  scale_null_shape_dist <- -2*(scale_loglik_from_scale - shape_loglik_from_scale)
  
  #shape(null)
  shape_null_constant_dist <- -2*(shape_loglik_from_shape - constant_loglik_from_shape)
  shape_null_scale_shape_dist <- -2*(shape_loglik_from_shape - scale_shape_loglik_from_shape)
  shape_null_scale_dist <- -2*(shape_loglik_from_shape - scale_loglik_from_shape)
  
  
  ##test statistic
  cat("Calculating Test Statistic. \n")
  constant_null_scale_shape_test <- 2*(gpd_constant$result$value - best_model_scale_shape$result$value)
  constant_null_scale_test <- 2*(gpd_constant$result$value - best_model_scale$result$value)
  constant_null_shape_test <- 2*(gpd_constant$result$value - best_model_shape$result$value)
  
  scale_shape_null_constant_test <- 2*(best_model_scale_shape$result$value - gpd_constant$result$value)
  scale_shape_null_scale_test <- 2*(best_model_scale_shape$result$value - best_model_scale$result$value)
  scale_shape_null_shape_test <- 2*(best_model_scale_shape$result$value - best_model_shape$result$value)
  
  scale_null_constant_test <- 2*(best_model_scale$result$value - gpd_constant$result$value)
  scale_null_scale_shape_test <- 2*(best_model_scale$result$value - best_model_scale_shape$result$value)
  scale_null_shape_test <- 2*(best_model_scale$result$value - best_model_shape$result$value)
  
  shape_null_constant_test <- 2*(best_model_shape$result$value - gpd_constant$result$value)
  shape_null_scale_shape_test <- 2*(best_model_shape$result$value - best_model_scale_shape$result$value)
  shape_null_scale_test <- 2*(best_model_shape$result$value - best_model_scale$result$value)
  
  ##Empirical P Value
  cat("Calculating P Value. \n")
  constant_null_scale_shape_pvalue <- sum(constant_null_scale_shape_dist > constant_null_scale_shape_test)/n_bootstrap
  constant_null_scale_pvalue <- sum(constant_null_scale_dist > constant_null_scale_test)/n_bootstrap
  constant_null_shape_pvalue <- sum(constant_null_shape_dist > constant_null_shape_test)/n_bootstrap
  
  scale_shape_null_constant_pvalue <- sum(scale_shape_null_constant_dist > scale_shape_null_constant_test)/n_bootstrap
  scale_shape_null_scale_pvalue <- sum(scale_shape_null_scale_dist > scale_shape_null_scale_test)/n_bootstrap
  scale_shape_null_shape_pvalue <- sum(scale_shape_null_shape_dist > scale_shape_null_shape_test)/n_bootstrap
  
  scale_null_constant_pvalue <- sum(scale_null_constant_dist > scale_null_constant_test)/n_bootstrap
  scale_null_scale_shape_pvalue <- sum(scale_null_scale_shape_dist > scale_null_scale_shape_test)/n_bootstrap
  scale_null_shape_pvalue <- sum(scale_null_shape_dist > scale_null_shape_test)/n_bootstrap
  
  shape_null_constant_pvalue <- sum(shape_null_constant_dist > shape_null_constant_test)/n_bootstrap
  shape_null_scale_shape_pvalue <- sum(shape_null_scale_shape_dist > shape_null_scale_shape_test)/n_bootstrap
  shape_null_scale_pvalue <- sum(shape_null_scale_dist > shape_null_scale_test)/n_bootstrap
  
  #
  pvalue_mat <- matrix(NA, 4, 4, dimnames = list(c("Constant(NULL)", "Scale_Shape(NULL)","Scale(NULL)", "Shape(NULL)"), c("Constant", "Scale_Shape","Scale", "Shape")))
  pvalue_mat[1,2] <- constant_null_scale_shape_pvalue
  pvalue_mat[1,3] <- constant_null_scale_pvalue
  pvalue_mat[1,4] <- constant_null_shape_pvalue
  
  pvalue_mat[2,1] <- scale_shape_null_constant_pvalue
  pvalue_mat[2,3] <- scale_shape_null_scale_pvalue
  pvalue_mat[2,4] <- scale_shape_null_shape_pvalue
  
  pvalue_mat[3,1] <- scale_null_constant_pvalue
  pvalue_mat[3,2] <- scale_null_scale_shape_pvalue
  pvalue_mat[3,4] <- scale_null_shape_pvalue
  
  pvalue_mat[4,1] <- shape_null_constant_pvalue
  pvalue_mat[4,2] <- shape_null_scale_shape_pvalue
  pvalue_mat[4,3] <- shape_null_scale_pvalue
                       
  write.csv(pvalue_mat, paste0(result_loc, "/Pvalue_mat.csv"), row.names = F)
  ##plot
  if(dir.exists(paste0(result_loc,"/constant_null"))==F) dir.create(paste0(result_loc,"/constant_null"))
  if(dir.exists(paste0(result_loc,"/scale_shape_null"))==F) dir.create(paste0(result_loc,"/scale_shape_null"))
  if(dir.exists(paste0(result_loc,"/scale_null"))==F) dir.create(paste0(result_loc,"/scale_null"))
  if(dir.exists(paste0(result_loc,"/shape_null"))==F) dir.create(paste0(result_loc,"/shape_null"))
  
  
  cat("Saving Plots(histogram of p value). \n")
  #constant null
  png(paste0(result_loc, "/constant_null/scale_shape.png"), 500, 500)
  hist(constant_null_scale_shape_dist, main = "Scale, Shape", 
       breaks=15,
       xlim = c(min(c(constant_null_scale_shape_test, constant_null_scale_shape_dist)), max(c(constant_null_scale_shape_test, constant_null_scale_shape_dist))))
  abline(v=constant_null_scale_shape_test, col = 'red')
  dev.off()
  
  png(paste0(result_loc, "/constant_null/scale.png"), 500, 500)
  hist(constant_null_scale_dist, main = "Scale",
       breaks=15,
       xlim = c(min(c(constant_null_scale_test, constant_null_scale_dist)), max(c(constant_null_scale_test, constant_null_scale_dist))))
  abline(v=constant_null_scale_test, col = 'red')
  dev.off()
  
  png(paste0(result_loc, "/constant_null/shape.png"), 500, 500)
  hist(constant_null_shape_dist, main = "Shape",
       breaks=15,
       xlim = c(min(c(constant_null_shape_test, constant_null_shape_dist)), max(c(constant_null_shape_test, constant_null_shape_dist))))
  abline(v=constant_null_shape_test, col = 'red')
  dev.off()
  
  #scale shape null
  png(paste0(result_loc, "/scale_shape_null/constant.png"), 500, 500)
  hist(scale_shape_null_constant_dist, main = "Constant",
       breaks=15,
       xlim = c(min(c(scale_shape_null_constant_test, scale_shape_null_constant_dist)), max(c(scale_shape_null_constant_test, scale_shape_null_constant_dist))))
  abline(v=scale_shape_null_constant_test, col = 'red')
  dev.off()
  
  png(paste0(result_loc, "/scale_shape_null/scale.png"), 500, 500)
  hist(scale_shape_null_scale_dist, main = "Scale",
       breaks=15,
       xlim = c(min(c(scale_shape_null_scale_test, scale_shape_null_scale_dist)), max(c(scale_shape_null_scale_test, scale_shape_null_scale_dist))))
  abline(v=scale_shape_null_scale_test, col = 'red')
  dev.off()
  
  png(paste0(result_loc, "/scale_shape_null/shape.png"), 500, 500)
  hist(scale_shape_null_shape_dist, main = "Shape",
       breaks=15,
       xlim = c(min(c(scale_shape_null_shape_test, scale_shape_null_shape_dist)), max(c(scale_shape_null_shape_test, scale_shape_null_shape_dist))))
  abline(v=scale_shape_null_shape_test, col = 'red')
  dev.off()
  
  #scale null
  png(paste0(result_loc, "/scale_null/constant.png"), 500, 500)
  hist(scale_null_constant_dist, main = "Constant",
       breaks=15,
       xlim = c(min(c(scale_null_constant_test, scale_null_constant_dist)), max(c(scale_null_constant_test, scale_null_constant_dist))))
  abline(v=scale_null_constant_test, col = 'red')
  dev.off()
  
  png(paste0(result_loc, "/scale_null/scale_shape.png"), 500, 500)
  hist(scale_null_scale_shape_dist, main = "Scale, Shape",
       breaks=15,
       xlim = c(min(c(scale_null_scale_shape_test, scale_null_scale_shape_dist)), max(c(scale_null_scale_shape_test, scale_null_scale_shape_dist))))
  abline(v=scale_null_scale_shape_test, col = 'red')
  dev.off()
  
  png(paste0(result_loc, "/scale_null/shape.png"), 500, 500)
  hist(scale_null_shape_dist, main = "Shape",
       breaks=15,
       xlim = c(min(c(scale_null_shape_test, scale_null_shape_dist)), max(c(scale_null_shape_test, scale_null_shape_dist))))
  abline(v=scale_null_shape_test, col = 'red')
  dev.off()
  
  #shape null
  png(paste0(result_loc, "/shape_null/constant.png"), 500, 500)
  hist(shape_null_constant_dist, main = "Constant",
       breaks=15,
       xlim = c(min(c(shape_null_constant_test, shape_null_constant_dist)), max(c(shape_null_constant_test, shape_null_constant_dist))))
  abline(v=shape_null_constant_test, col = 'red')
  dev.off()
  
  png(paste0(result_loc, "/shape_null/scale_shape.png"), 500, 500)
  hist(shape_null_scale_shape_dist, main = "Scale, Shape",
       breaks=15,
       xlim = c(min(c(shape_null_scale_shape_test, shape_null_scale_shape_dist)), max(c(shape_null_scale_shape_test, shape_null_scale_shape_dist))))
  abline(v=constant_null_scale_shape_test, col = 'red')
  dev.off()
  
  png(paste0(result_loc, "/shape_null/scale.png"), 500, 500)
  hist(shape_null_scale_dist, main = "scale",
       breaks=15,
       xlim = c(min(c(shape_null_scale_test, shape_null_scale_dist)), max(c(shape_null_scale_test, shape_null_scale_dist))))
  abline(v=constant_null_scale_test, col = 'red')
  dev.off()
  
  ##RETURN result
  cat("Finishing Bootstrap\n")
  result_list <- list(mle_mat_scale_shape = mle_mat_scale_shape,
                      mle_mat_scale = mle_mat_scale,
                      mle_mat_shape = mle_mat_shape,
                      
                      gpd_constant = gpd_constant,
                      best_model_scale_shape = best_model_scale_shape,
                      best_model_scale = best_model_scale,
                      best_model_shape = best_model_shape,
                      #Log likelihood
                      constant_loglik_from_constant = constant_loglik_from_constant,
                      constant_loglik_from_scale_shape = constant_loglik_from_scale_shape,
                      constant_loglik_from_scale = constant_loglik_from_scale,
                      constant_loglik_from_shape = constant_loglik_from_shape,
                      
                      scale_shape_loglik_from_constant = scale_shape_loglik_from_constant,
                      scale_shape_loglik_from_scale_shape = scale_shape_loglik_from_scale_shape,
                      scale_shape_loglik_from_scale = scale_shape_loglik_from_scale,
                      scale_shape_loglik_from_shape = scale_shape_loglik_from_shape,
                      
                      scale_loglik_from_constant = scale_loglik_from_constant,
                      scale_loglik_from_scale_shape = scale_loglik_from_scale_shape,
                      scale_loglik_from_scale = scale_loglik_from_scale,
                      scale_loglik_from_shap = scale_loglik_from_shape,
                      
                      shape_loglik_from_constant = shape_loglik_from_constant,
                      shape_loglik_from_scale_shape = shape_loglik_from_scale_shape,
                      shape_loglik_from_scale = shape_loglik_from_scale,
                      shape_loglik_from_shape = shape_loglik_from_shape,
                      
                      #TEST STATISTIC
                      constant_null_scale_shape_test = constant_null_scale_shape_test,
                      constant_null_scale_test = constant_null_scale_test,
                      constant_null_shape_test = constant_null_shape_test,
                      
                      scale_shape_null_constant_test = scale_shape_null_constant_test,
                      scale_shape_null_scale_test = scale_shape_null_scale_test,
                      scale_shape_null_shape_test = scale_shape_null_shape_test,
                      
                      scale_null_constant_test = scale_null_constant_test,
                      scale_null_scale_shape_test = scale_null_scale_shape_test,
                      scale_null_shape_test = scale_null_shape_test,
                      
                      shape_null_constant_test = shape_null_constant_test,
                      shape_null_scale_shape_test = shape_null_scale_shape_test,
                      shape_null_scale_test = shape_null_scale_test,
                      #P VALUE
                      constant_null_scale_shape_pvalue = constant_null_scale_shape_pvalue,
                      constant_null_scale_pvalue = constant_null_scale_pvalue,
                      constant_null_shape_pvalue = constant_null_shape_pvalue,
                      
                      scale_shape_null_constant_pvalue = scale_shape_null_constant_pvalue,
                      scale_shape_null_scale_pvalue = scale_shape_null_scale_pvalue,
                      scale_shape_null_shape_pvalue = scale_shape_null_shape_pvalue,
                      
                      scale_null_constant_pvalue = scale_null_constant_pvalue,
                      scale_null_scale_shape_pvalue = scale_null_scale_shape_pvalue,
                      scale_null_shape_pvalue = scale_null_shape_pvalue,
                      
                      shape_null_constant_pvalue = shape_null_constant_pvalue,
                      shape_null_scale_shape_pvalue = shape_null_scale_shape_pvalue,
                      shape_null_scale_pvalue = shape_null_scale_pvalue
  )
  save(result_list, file= paste0(result_loc, "/result_list.RData"))
  return(result_list)
}
