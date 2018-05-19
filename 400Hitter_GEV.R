# install.packages("Lahman")
# install.packages("dplyr")
# install.packages("ggplot2")
# install.packages("extRemes")
install.packages(c("Lahman", "dplyr", "ggplot2", "extRemes","splines2", "fExtremes"))
# install.packages("ismev")
# install.packages("splines2")
library(splines2)
library(ismev)
library(ggplot2)
library(dplyr)
library(Lahman)
library(extRemes)
library(fExtremes)
library(splines)

#########################################################################################
#################################DATA HANDLING###########################################
#########################################################################################
#Load data
# data(package = 'Lahman')
data("Batting")

#Data Description
#G : Games Played
#AB : Bats
#R : Run
#H : Hits()
# X2B : Doulbs hits on which the batter reached second base safely
# X3B : Triples: hits on which the batter reached third base safely
# HR : Homeruns
# RBI : Runs Batted In
# SB : Stolen Bases
# CS : Caught Stealing
# BB : Base on Balls
# SO : Strikeouts
# IBB : Intentional walks
# HBP : Hit by pitch
# SH : Sacrifice Hits
# SF : Sacrifice files
# GIDP : Grounded into double plays

Lahman::Salaries
#Data Handling
batting <- battingStats()
#add salary
batting <- merge(batting, 
                 Lahman::Salaries[,c("playerID", "yearID", "teamID", "salary")], 
                 by=c("playerID", "yearID", "teamID"), all.x=TRUE)

#Add name, age and bat hand information:
masterInfo <- Lahman::Master[, c('playerID', 'birthYear', 'birthMonth',
                                 'nameLast', 'nameFirst', 'bats')]
batting <- merge(batting, masterInfo, all.x = TRUE)
batting$age <- with(batting, yearID - birthYear -
                      ifelse(birthMonth < 10, 0, 1))
#sorting
batting <- arrange(batting, playerID, yearID, stint)

#filter for eligible hitters over PA > 450 after 1900.
eligibleHitters <- batting %>% dplyr::filter(yearID >= 1900 & PA > 450)

topHitters  <- eligibleHitters %>% group_by(yearID) %>%dplyr::filter(BA > .350)

# Create a factor variable to distinguish the .400 hitters
topHitters$ba400 <- with(topHitters, BA >= 0.400)

# Sub-data frame for the .400 hitters plus the outliers after 1950
# (averages above .380) - used to produce labels in the plot below
bignames <- rbind(subset(topHitters, ba400),
                  subset(topHitters, yearID > 1950 & BA > 0.380))
# Cut to the relevant set of variables
bignames <- subset(bignames, select = c('playerID', 'yearID', 'nameLast',
                                        'nameFirst', 'BA'))

# Positional offsets to spread out certain labels
#                     NL TC JJ TC GS TC RH GS HH RH RH BT TW TW  RC GB TG 
bignames$xoffset <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, -2, 0, 0)
bignames$yoffset <- c(0, 0, -0.003, 0, 0, 0, 0, 0, -0.004, 0, 0, 0, 0, 0, -0.003, 0, 0)  +  0.002



ggplot(topHitters, aes(x = yearID, y = BA)) +
  geom_point(aes(colour = ba400), size = 2.5) +
  geom_hline(yintercept = 0.400, size = 1) +
  geom_text(data = bignames, aes(x = yearID + xoffset, y = BA + yoffset,
                                 label = nameLast), size = 3) +
  ylim(0.330, 0.430) +
  xlab('Year') +
  scale_y_continuous('Batting average',
                     breaks = seq(0.34, 0.42, by = 0.02),
                     labels = c('.340', '.360', '.380', '.400', '.420')) +
  geom_smooth(method = "loess")


#Pitcher
Lahman::Pitching
#Data Handling
batting <- battingStats()
#add salary
batting <- merge(batting, 
                 Lahman::Salaries[,c("playerID", "yearID", "teamID", "salary")], 
                 by=c("playerID", "yearID", "teamID"), all.x=TRUE)

#Add name, age and bat hand information:
masterInfo <- Lahman::Master[, c('playerID', 'birthYear', 'birthMonth',
                                 'nameLast', 'nameFirst', 'bats')]
batting <- merge(batting, masterInfo, all.x = TRUE)
batting$age <- with(batting, yearID - birthYear -
                      ifelse(birthMonth < 10, 0, 1))
#sorting
batting <- arrange(batting, playerID, yearID, stint)

#use eligibleHitters$BA as data in GPD model
GPD_data_temp <- eligibleHitters[!is.na(eligibleHitters$BA),]

#########################################################################################
#############################Fitting GPD with extRemes###################################
#########################################################################################
## GEV fitting
gev_fit_model <- extRemes::fevd(topHitters$BA, type = "GEV", method = "MLE")
summary(gev_fit_model)
plot(gev_fit_model)


GPD_data_temp$BA
####
####Threshold
threshrange.plot(GPD_data_temp$BA, r = c(0.3, 0.4), nint = 11)

#constant threshold
thrhld <- 0.34
GPD_data <- GPD_data_temp
GPD_data <- GPD_data_temp[GPD_data_temp$BA > thrhld,]
plot(GPD_data$year, GPD_data$BA)

#generate non-constant threshold
top_n_threshold <- aggregate(BA~yearID, GPD_data_temp, function(x) sort(x,decreasing =T)[5+1])

threshold_of_year<-rep(NA, nrow(GPD_data_temp))

for(i in 1:nrow(GPD_data_temp)){
  threshold_of_year[i] <- top_n_threshold[which(GPD_data_temp[i,"yearID"] == top_n_threshold$yearID),"BA"]
}

GPD_data <- GPD_data_temp[GPD_data_temp$BA > threshold_of_year,]


#GPD mean function
GPD_mean <- function(mu, sigma, xi){
  return(mu + sigma/(1-xi))
}


GPD_data <- GPD_data_temp[GPD_data_temp$BA > threshold,]

iSpline_list <- vector("list", 16)
mle_mat <- matrix(NA, ncol = 4, nrow = 4)

x_range <- GPD_data$year[GPD_data$BA > thrhld]
threshold <- 0.37
count <- 0
for(i in 1:4){
  for(j in 1:4){
    cat("i = ", i, ", j = ", j, "\n")
    count <- count + 1
    n_knots <- i
    dgr <- j
    if(i==1){
      knots <- round(quantile(GPD_data$year,0.5))
    }else if(i==2){
      knots <- c(round(quantile(GPD_data$year,0.33)), round(quantile(GPD_data$year,0.66)))
    }else if(i==3){
      knots <- c(round(quantile(GPD_data$year,0.25)), round(quantile(GPD_data$year,0.5)), round(quantile(GPD_data$year,0.75)))
    }else if(i==4){
      knots <- c(round(quantile(GPD_data$year,0.2)), round(quantile(GPD_data$year,0.4)), round(quantile(GPD_data$year,0.6)), round(quantile(GPD_data$year,0.8)))
    }
    
    #fitting
    gpd_ispline <- fevd(GPD_data$BA, threshold = thrhld, type = "GP", method = "MLE",
                        shape.fun = ~ iSpline(c(1899,GPD_data$year),knots=knots,degree=dgr,intercept=T)[-1,]-1,
                        scale.fun = ~ iSpline(c(1899,GPD_data$year),knots=knots,degree=dgr,intercept=T)[-1,]-1
    )
    
    iSpline_list[[count]] <- gpd_ispline
    names(iSpline_list[count]) <- paste0("knots = ", i, ", degree = ", j, "\n")
    
    mle_mat[i,j] <- gpd_ispline$results$value
    
    # iSpline_basis <- iSpline(c(1899,GPD_data$year),knots=knots,degree=dgr,intercept=T)[-1,]
    
    #obtain parameters for spline
    # sigma_par <- gpd_ispline$results$par[1:(length(gpd_ispline$results$par)-1)]
    # sigma_par <- gpd_ispline$results$par[1:(length(gpd_ispline$results$par)/2)]
    # xi_par <- gpd_ispline$results$par[((length(gpd_ispline$results$par)/2) + 1):length(gpd_ispline$results$par)]
    
    #predict GPD parameter vector
    # sigma_vec <- iSpline_basis%*%as.matrix(sigma_par)
    # xi_vec <- iSpline_basis%*%as.matrix(xi_par)
    # xi_vec <- rep(gpd_ispline$results$par[length(gpd_ispline$results$par)], length(sigma_vec))
    
    sigma_vec <- findpars(gpd_ispline)$scale
    xi_vec <- findpars(gpd_ispline)$shape
    GPD_mean_vec <- GPD_mean(mu = threshold, sigma = sigma_vec, xi = xi_vec)
    
    #plot
    png(paste0("C:/Users/seho1/Documents/parameter_plot/",i,"knots, ", j, "degree.png"), 500, 500)
    par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
    plot(x_range, sigma_vec, xlab = "year", ylab = "sigma", main = "year - sigma")
    plot(x_range, xi_vec, xlab = "year", ylab = "xi", main = "year-xi")
    plot(x_range, GPD_mean_vec, xlab = "year", ylab = "GPD mean", main = "year - GPD mean")
    # mtext(paste0(), outer = TRUE, cex = 1.5)
    dev.off()
    par(mfrow=c(1,1))
  }
}

# write.csv(mle_mat, "mle_mat.csv", row.names = F)


###############################################################
#sample from constant model
degree=4
gpd_constant <- fevd(GPD_data$BA, threshold = thrhld, type = "GP", method = "MLE")
# knots <- c(round(quantile(GPD_data$year,0.33)), round(quantile(GPD_data$year,0.66)))
# knots <- c(round(quantile(GPD_data$year,0.25)), round(quantile(GPD_data$year,0.5)), round(quantile(GPD_data$year,0.75)))
knots <- c(round(quantile(GPD_data$year,0.2)), round(quantile(GPD_data$year,0.4)), round(quantile(GPD_data$year,0.6)), round(quantile(GPD_data$year,0.8)))
gpd_ispline <- fevd(GPD_data$BA, threshold = thrhld, type = "GP", method = "MLE",
                    shape.fun = ~ iSpline(c(1899,GPD_data$year),knots=knots,degree=4,intercept=T)[-1,]-1,
                    scale.fun = ~ iSpline(c(1899,GPD_data$year),knots=knots,degree=4,intercept=T)[-1,]-1)




# gpd_ispline <- iSpline_list[[14]]
# iSpline_basis <- eval(parse(text=paste0("iSpline(GPD_data$year,knots=knots,degree=",4,",intercept=T)")))[GPD_data$BA > 0.37,]
# iSpline_basis <- iSpline(c(1899,GPD_data$year),knots=knots,degree=4,intercept=T)[-1,]
#obtain parameters for spline
# sigma_par <- gpd_ispline$results$par[1:(length(gpd_ispline$results$par)-1)]
# sigma_par <- gpd_ispline$results$par[1:(length(gpd_ispline$results$par)/2)]
# xi_par <- gpd_ispline$results$par[((length(gpd_ispline$results$par)/2) + 1):length(gpd_ispline$results$par)]
#predict GPD parameter vector
# xi_vec <- rep(as.numeric(gpd_ispline$results$par[length(gpd_ispline$results$par)]), length(sigma_vec))

sigma_vec <- findpars(gpd_ispline)$scale
xi_vec <- findpars(gpd_ispline)$shape

n_bootstrap <- 1000
null_loglik <- rep(NA, n_bootstrap)
alter_loglik <- rep(NA, n_bootstrap)
for(i in 1:n_bootstrap){
  cat(i, "\n")
  bootstrap_data <- data.frame(year = GPD_data$yearID[GPD_data$BA>thrhld], BA = rextRemes(x=gpd_constant, sum(GPD_data$BA>thrhld)))
  #null
  null_gpd_constant <- fevd(bootstrap_data2$BA, threshold = thrhld, type = "GP", method = "MLE")
  # null_loglik[i] <- levd(bootstrap_data[,2], threshold = thrhld, scale = gpd_constant$results$par[1], shape = gpd_constant$results$par[2], type = "GP", negative=F)
  # log(prod(d_gpd(x= bootstrap_data[,2], sigma = gpd_constant$results$par[1], xi = gpd_constant$results$par[2], threshold = 0.37)))
  null_loglik[i] <- -null_gpd_constant$results$value
  #alternative
  alter_gpd_ispline <- fevd(bootstrap_data2$BA, threshold = thrhld, type = "GP", method = "MLE",
                            shape.fun = ~ iSpline(c(1899,GPD_data$year),knots=knots,degree=4,intercept=T)[-1,]-1,
                            scale.fun = ~ iSpline(c(1899,GPD_data$year),knots=knots,degree=4,intercept=T)[-1,]-1)
  alter_loglik[i] <- -alter_gpd_ispline$result$value
  # alter_loglik[i] <- levd(bootstrap_data[,2], threshold = 0.37, scale = as.numeric(sigma_vec), shape = as.numeric(xi_vec), type = "GP", negative = F)
  # alter_loglik[i] <- levd(bootstrap_data[,2], threshold = thrhld, scale = as.numeric(sigma_vec), shape = as.numeric(xi_vec), type = "GP", negative = F)
}


#test statistics
null_test_stat <- 2*(alter_loglik - null_loglik)
test_stat <- 2*(-gpd_ispline$results$value + gpd_constant$results$value)
empirical_pvalue <- sum(null_test_stat >test_stat)/1000


####iSpline(null) - constant(alternative)
n_bootstrap <- 1000
null_loglik2 <- rep(NA, n_bootstrap)
alter_loglik2 <- rep(NA, n_bootstrap)

gpd_ispline_pars <- findpars(gpd_ispline)

for(i in 1:n_bootstrap){
  cat(i, "\n")
  BA <- rep(NA, nrow(GPD_data))
  for(j in 1:nrow(GPD_data)){
    BA[j] <- revd(n=1, scale=gpd_ispline_pars$scale[j], shape=gpd_ispline_pars$shape[j], threshold =thrhld, type="GP")
  }
  bootstrap_data2 <- data.frame(year = GPD_data$yearID[GPD_data$BA>thrhld], BA = BA)
  #null
  null_gpd_ispline <- fevd(bootstrap_data$BA, threshold = thrhld, type = "GP", method = "MLE",
                           shape.fun = ~ iSpline(c(1899,GPD_data$year),knots=knots,degree=4,intercept=T)[-1,]-1,
                           scale.fun = ~ iSpline(c(1899,GPD_data$year),knots=knots,degree=4,intercept=T)[-1,]-1)
  null_loglik2[i] <- -null_gpd_ispline$result$value
  
  #alternative
  alter_gpd_constant <- fevd(bootstrap_data$BA, threshold = thrhld, type = "GP", method = "MLE")
  # null_loglik[i] <- levd(bootstrap_data[,2], threshold = thrhld, scale = gpd_constant$results$par[1], shape = gpd_constant$results$par[2], type = "GP", negative=F)
  # log(prod(d_gpd(x= bootstrap_data[,2], sigma = gpd_constant$results$par[1], xi = gpd_constant$results$par[2], threshold = 0.37)))
  alter_loglik2[i] <- -alter_gpd_constant$results$value
  
  # alter_loglik[i] <- levd(bootstrap_data[,2], threshold = 0.37, scale = as.numeric(sigma_vec), shape = as.numeric(xi_vec), type = "GP", negative = F)
  # alter_loglik[i] <- levd(bootstrap_data[,2], threshold = thrhld, scale = as.numeric(sigma_vec), shape = as.numeric(xi_vec), type = "GP", negative = F)
}


#test statistics
null_test_stat <- 2*(alter_loglik2 - null_loglik2)
test_stat <- 2*(+gpd_ispline$results$value - gpd_constant$results$value)
empirical_pvalue2 <- sum(null_test_stat >test_stat)/1000

a <- d_gpd(bootstrap_data[,2], sigma = as.numeric(sigma_vec), xi = as.numeric(xi_vec), threshold=thrhld)

d_gpd<- function(x, sigma, xi, threshold){
  return((1/sigma)*(1+(xi*(x-threshold))/sigma)^(-1/xi - 1))
}


#####################################################################################
#####################################################################################
#Natural Cubic Spline
gpd_ncspline <- fevd(GPD_data$BA, threshold = 0.37, type = "GP", method = "MLE",
                     scale.fun = ~ ns(GPD_data$year,df = 2,intercept=F),
                     shape.fun = ~ ns(GPD_data$year,df = 2,intercept=F)
)

ns(GPD_data$year,df = dof,intercept=T)[order(GPD_data$year),]

n <- 10
ncs_list <- vector("list", n)
ncs_mle_mat <- matrix(NA, ncol = 1, nrow = n)

x_range <- GPD_data$year[GPD_data$BA > 0.37]
threshold <- 0.37
count <- 0
for(i in 1:n){
  cat("i = ", i,"\n")
  count <- count + 1
  dof <- i
  
  #fitting
  gpd_ncspline <- fevd(GPD_data$BA, threshold = 0.37, type = "GP", method = "MLE",
                       shape.fun = ~ ns(GPD_data$year,df = dof,intercept=F)-1,
                       scale.fun = ~ ns(GPD_data$year,df = dof,intercept=F)-1
  )
  
  ncs_list[[count]] <- gpd_ncspline
  names(ncs_list[count]) <- paste0("degree of freedome = ", i, "\n")
  
  ncs_mle_mat[i,1] <- gpd_ncspline$results$value
  
  ncspline_basis <- ns(GPD_data$year,df = dof,intercept=T)
  
  #obtain parameters for spline
  # sigma_par <- gpd_ispline$results$par[1:(length(gpd_ispline$results$par)-1)]
  sigma_par <- gpd_ncspline$results$par[1:(length(gpd_ncspline$results$par)/2)]
  xi_par <- gpd_ncspline$results$par[((length(gpd_ncspline$results$par)/2) + 1):length(gpd_ncspline$results$par)]
  
  #predict GPD parameter vector
  sigma_vec <- ncspline_basis%*%as.matrix(sigma_par)
  xi_vec <- ncspline_basis%*%as.matrix(xi_par)
  # xi_vec <- rep(gpd_ispline$results$par[length(gpd_ispline$results$par)], length(sigma_vec))
  GPD_mean_vec <- GPD_mean(mu = threshold, sigma = sigma_vec, xi = xi_vec)
  
  #plot
  png(paste0("C:/Users/seho1/Documents/parameter_plot/ncs/",i,"df.png"), 500, 500)
  par(mfrow=c(2,2),oma = c(0, 0, 2, 0))
  plot(x_range, sigma_vec, xlab = "year", ylab = "sigma", main = "year - sigma")
  plot(x_range, xi_vec, xlab = "year", ylab = "xi", main = "year-xi")
  plot(x_range, GPD_mean_vec, xlab = "year", ylab = "GPD mean", main = "year - GPD mean")
  # mtext(paste0(), outer = TRUE, cex = 1.5)
  dev.off()
  par(mfrow=c(1,1))
  
}
iSpline_list[[2]]
findpars(iSpline_list[[14]])



ncol(GPD_data) - (1900 - 2016 + 1) #151
knots <- c(round(quantile(GPD_data$year,0.25)), round(quantile(GPD_data$year,0.5)), round(quantile(GPD_data$year,0.75)))
gpd_const <- fevd(GPD_data$BA, threshold = 0.37, type = "GP")
gpd_ispline <- fevd(GPD_data$BA, threshold = 0.37, type = "GP", method = "MLE",
                    scale.fun = ~ iSpline(GPD_data$year,knots=knots,degree=2,intercept=F),
                    shape.fun = ~ iSpline(GPD_data$year,knots=knots,degree=2,intercept=F)
)
gpd_ncspline <- fevd(GPD_data$BA, threshold = 0.37, type = "GP", method = "MLE",
                     scale.fun = ~ ns(GPD_data$year,df = 2,intercept=F),
                     shape.fun = ~ ns(GPD_data$year,df = 2,intercept=F)
)

#ncs
NCS_basis<- ns(GPD_data$year[GPD_data$BA >= 0.37],df = 3,intercept=F)
sigma_par <- gpd_ncspline$results$par[1:(length(gpd_ncspline$results$par)/2)]
xi_par <- gpd_ncspline$results$par[((length(gpd_ncspline$results$par)/2) + 1):length(gpd_ncspline$results$par)]
a <- NCS_basis%*%as.matrix(sigma_par)
plot(GPD_data$year[GPD_data$BA >= 0.37], a)
