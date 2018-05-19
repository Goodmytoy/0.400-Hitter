GPD_data_temp <- eligibleHitters[!is.na(eligibleHitters$BA),]

#iter = 100
gpdtest_0.34_int_100 <- GPD_sim_test(GPD_data_temp, threshold = 0.34,intercept=T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.34_int_100", n_bootstrap = 100)
gpdtest_0.35_int_100 <- GPD_sim_test(GPD_data_temp, threshold = 0.35,intercept= T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.35_int_100", n_bootstrap = 100, iter_unit=100)
gpdtest_0.36_int_100 <- GPD_sim_test(GPD_data_temp, threshold = 0.36,intercept= T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.36_int_100", n_bootstrap = 100)
gpdtest_int_top5_100 <- GPD_sim_test(GPD_data_temp, top_n = 5, intercept=T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_int_top5_100", n_bootstrap = 100)
#iter = 1000
gpdtest_0.34_int_1000 <- GPD_sim_test(GPD_data_temp, threshold = 0.34,intercept=T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.34_int_1000", n_bootstrap = 1000)
gpdtest_0.35_int_1000 <- GPD_sim_test(GPD_data_temp, threshold = 0.35,intercept= T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.35_int_1000", n_bootstrap = 1000, iter_unit=100)
gpdtest_0.36_int_1000 <- GPD_sim_test(GPD_data_temp, threshold = 0.36,intercept= T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.36_int_1000", n_bootstrap = 1000)
gpdtest_int_top5_1000 <- GPD_sim_test(GPD_data_temp, top_n = 5, intercept=T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_int_top5_1000", n_bootstrap = 1000)
#iter = 10000
gpdtest_0.34_int_10000 <- GPD_sim_test(GPD_data_temp, threshold = 0.34,intercept=T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.34_int_10000", n_bootstrap = 10000)
gpdtest_0.35_int_10000 <- GPD_sim_test(GPD_data_temp, threshold = 0.35,intercept= T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.35_int_10000", n_bootstrap = 10000, iter_unit=100)
gpdtest_0.36_int_10000 <- GPD_sim_test(GPD_data_temp, threshold = 0.36,intercept= T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.36_int_10000", n_bootstrap = 10000)
gpdtest_int_top5_10000 <- GPD_sim_test(GPD_data_temp, top_n = 5, intercept=T, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_int_top5_10000", n_bootstrap = 10000)

# gpdtest_0.35 <- GPD_sim_test(GPD_data_temp, threshold = 0.35,intercept=F, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.35", n_bootstrap = 10000)
# gpdtest_0.34 <- GPD_sim_test(GPD_data_temp, threshold = 0.34,intercept=F, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.34", n_bootstrap = 10000)
# gpdtest_0.36 <- GPD_sim_test(GPD_data_temp, threshold = 0.36,intercept=F, result_loc = "D:/문서 백업/박세호&강상욱 논문/GPDtest_0.36", n_bootstrap = 10000)
const_loglik <- -gpdtest_0.34_int$gpd_constant$results$value
scale_shape_loglik <- -gpdtest_0.34_int$best_model_scale_shape$results$value

-2*(const_loglik - scale_shape_loglik)
-2*(scale_shape_loglik - const_loglik)

gpdtest_0.34_int_top5$constant_null_scale_shape_pvalue
gpdtest_0.34_int_top5$scale_shape_null_constant_pvalue


hist(-2*(gpdtest_0.34_int$constant_loglik_from_constant - gpdtest_0.34_int$scale_shape_loglik_from_constant))
abline(v=-2*(const_loglik - scale_shape_loglik), col = 2)
sum(gpdtest_0.34_int$constant_loglik_from_constant > gpdtest_0.34_int$scale_shape_loglik_from_constant)/10000

hist(-2*(gpdtest_0.34_int$scale_shape_loglik_from_scale_shape - gpdtest_0.34_int$constant_loglik_from_scale_shape))
abline(v=-2*(scale_shape_loglik - const_loglik), col = 2)

sum(gpdtest_0.34_int$constant_loglik_from_scale_shape > gpdtest_0.34_int$scale_shape_loglik_from_scale_shape)/10000

dev.off()



