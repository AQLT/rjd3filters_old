library(rjd3filters)
y = rjd3toolkit::ABS$X0.2.09.10.M
e1 <- moving_average(rep(1,12), lags = -6)
e1 <- e1/sum(e1)
e2 <- moving_average(rep(1/12, 12), lags = -5)
# used to have the 1rst estimate of the trend
tc_1 <- M2X12 <- (e1 + e2)/2
coef(M2X12) |> round(3)
si_1 <- 1 - tc_1

M3X3 <- rjd3filters::macurves("S3x3")
M3X3_s <- rjd3filters:::finite_filters(moving_average(M3X3[,3],-2),
                                      list(moving_average(M3X3[,1], -2),
                                           moving_average(M3X3[,2], -2)),first_to_last = T)
M3X3_s = to_seasonal(M3X3_s, 12)
as.matrix(M3X3_s) |>  rownames()
s_1 = M3X3_s * si_1
s_2 = M3X3_s * finite_filters(si_1)
e1 = M3X3_s
e2 = finite_filters(si_1)
mm_s_norm = finite_filters((1 - M2X12), rfilters = lapply(1:6, function(i) (1 - M2X12)*moving_average(1, lags = -i)))
e1 = mm_s_norm
e2 = s_1
s_1_norm = mm_s_norm * s_2
sa_1 <- 1 - s_1_norm
sa_1*y

henderson_mm = finite_filters(lp_filter(horizon = 6,ic = 3.5))
tc_2 <- henderson_mm * sa_1
si_2 <- 1 - tc_2
M5 <- moving_average(rep(1/5, 5), lags = -2)
M5X5_seasonal <- to_seasonal(M5 * M5, 12)
s_2 = M5X5_seasonal * si_2
s_2_norm = mm_s_norm * s_2
sa_2 <- 1 - s_2_norm
tc_f <- henderson_mm * sa_2
tc_f * y
tmp1 = sa_2* y
tmp2 = x11(y = y, trend.coefs = lp_filter(horizon = 6,ic = 3.5)$filters.coef,extreme.lsig = 300, extreme.usig = 400, mul = FALSE)
tmp1 - tmp2$decomposition[,"sa"]
plot(tmp1)
lines(tmp2$decomposition[,"sa"], col = "red")

# IC-ratio varie ici
tmp = rjd3x13::x11(y, rjd3x13::spec_x11_default() |>
               rjd3x13::set_x11(mode = "Additive",
                                seasonal.filter = "S3X5",
                                henderson.filter = 13,
                                lsigma = 3000,
                                usigma = 4000))
tmp2$decomposition[,"sa"]
tmp$d11
