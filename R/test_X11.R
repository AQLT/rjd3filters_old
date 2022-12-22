# library(rjd3filters)
# y = rjd3toolkit::ABS$X0.2.09.10.M
#
# # wk <- RJDemetra::new_workspace()
# # RJDemetra::new_multiprocessing(wk, "sa1")
# # RJDemetra::add_sa_item(wk, "sa1", RJDemetra::x13(y, spec = "X11"), "X13")
# # RJDemetra::save_workspace(wk, "workspace.xml")
# wk <- RJDemetra::load_workspace("workspace.xml")
# RJDemetra::compute(wk)
# mod <- RJDemetra::get_jmodel(wk,progress_bar = FALSE)[[1]][[1]]
# extract <- function(id) {
#   RJDemetra::get_indicators(mod, sprintf("decomposition.%s", id))[[1]]
# }
# compare <- function(x, id){
#   res = ts.intersect(x, extract(id))
#   all.equal(res[,1], res[,2])
# }
# # extract <- function(x, id){
# #   ts(rjd3toolkit:::proc_vector(x$java, id),
# #      end = end(y), frequency = frequency(y))
# # }
# # x11_step = x11(y = y, trend.coefs = lp_filter(horizon = 6,ic = 3.5)$filters.coef,extreme.lsig = 300, extreme.usig = 400, mul = FALSE)
# waldo::compare(y, extract("b1"))
#
# e1 <- moving_average(rep(1,12), lags = -6)
# e1 <- e1/sum(e1)
# e2 <- moving_average(rep(1/12, 12), lags = -5)
# # used to have the 1rst estimate of the trend
# tc_1 <- M2X12 <- (e1 + e2)/2
# coef(M2X12) |> round(3)
# compare(M2X12 * y, "b2")
# si_1 <- 1 - tc_1
# compare(si_1 * y, "b3")
#
# M3X3 <- macurves("S3x3")
# M3X3_s <- to_seasonal(M3X3, 12)
# as.matrix(M3X3_s) |>  rownames()
# as.matrix(M3X3, TRUE, TRUE)
# s_1 = M3X3_s * si_1
# mm_s_norm = finite_filters((1 - M2X12),
#                            rfilters = lapply(1:6, function(i) 1 - M2X12*moving_average(1, lags = -i))
#                            )
# impute_6_last <- finite_filters(moving_average(c(rep(0, 6), 1), lags = -6),
#                                 lapply(1:6, function(i) moving_average(1, lags = -12)))
#
# s_non_norm = M3X3_s * (si_1 * y)
# s_non_norm[is.na(s_non_norm)] <- 0
# mm_s_norm * s_non_norm
# round((1 - M2X12) * s_non_norm, 2) - jfilter(s_non_norm, mm_s_norm, remove_missing = FALSE) |> round(2)
# mm_s_norm@sfilter
# s_1_y = jfilter(M3X3_s * (si_1 * y), mm_s_norm, remove_missing = FALSE)
# s_1_y[is.na(s_1_y)] <- 0
# impute_6_last * s_1_y
#
# jfilter(y = s_1_y, coefs = impute_6_last, remove_missing = FALSE)
# y = s_1_y
# impute_6_last = impute_6_last
# compare(s_1_y, "b4")
#
#
# waldo::compare(M3X3_s * (si_1 * y), s_2 * y)
# AQLTools::hc_stocks(ts.union(M3X3_s * (si_1 * y), s_2 * y))
#
# compare(x = (1 - M2X12) * (M3X3_s * (si_1 * y)), "b5")
#
# AQLTools::hc_stocks(ts.union((1 - M2X12) * (M3X3_s * (si_1 * y)),
#                              extract("b5")))
#
#
# mm_s_norm = finite_filters((1 - M2X12), rfilters = lapply(1:6, function(i) (1 - M2X12) * moving_average(1, lags = -12-i)
# ))
# round(as.matrix(mm_s_norm),3)
# s_1_norm = (1 - M2X12) * s_2
# (1 - M2X12) * s_2 * y - (1 - M2X12) * (s_2 * y)
# s_2_norm = mm_s_norm * s_2
# round(as.matrix(s_2),3) |> View()
# # finite_filters(moving_average(rep(1/3,3), -1),
# #                rfilters = moving_average(c(1/3,1/3,NA), -1))
# round(as.matrix(s_1_norm),3) |> View()
# round(as.matrix(s_2_norm),3) |> View()
#
# compare(s_1_norm * y, "b5")
# waldo::compare(s_1_norm * y, s_2_norm * y)
# s_2_norm * y
# s_2_norm@sfilter@upper_bound
# AQLTools::hc_stocks(ts.union(s_1_norm * y, extract("b5")))
# sa_1 <- 1 - s_1_norm
# sa_1*y
#
# henderson_mm = finite_filters(lp_filter(horizon = 6,ic = 3.5))
# tc_2 <- henderson_mm * sa_1
# tc_2*y
#
# rjd3x13::x13()
#
# si_2 <- 1 - tc_2
# si_2*y
#
# M3X5 <- rjd3filters::macurves("S3x5")
# M3X5_s = to_seasonal(M3X5, 12)
# s_2 = M3X5_s * si_2
# dim(as.matrix(si_2))
# dim(as.matrix(s_2)) #pair, pb
# e1 = M3X5_s
# e2 = si_2
# M3X5_s * y
# si_2*y
# s_2_norm = mm_s_norm * s_2
# sa_2 <- 1 - s_2_norm
# tc_f <- henderson_mm * sa_2
# tc_f * y
# tmp1 = sa_2* y
#
# tmp1 - tmp2$decomposition[,"sa"]
# plot(tmp1)
# lines(tmp2$decomposition[,"sa"], col = "red")
#
# # IC-ratio varie ici
# tmp = rjd3x13::x11(y, rjd3x13::spec_x11_default() |>
#                      rjd3x13::set_x11(mode = "Additive",
#                                       seasonal.filter = "S3X5",
#                                       henderson.filter = 13,
#                                       lsigma = 3000,
#                                       usigma = 4000))
# tmp2$decomposition[,"sa"]
# tmp$d11
#
