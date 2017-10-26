
# Options -----------------------------------------------------------------

library("SemiCompRisks")
library("scales")

set.seed(42)



# Data simulation ---------------------------------------------------------

# Data generation parameters
n <- 1500
beta1.true <- c(0.1, 0.1)
beta2.true <- c(0.2, 0.2)
beta3.true <- c(0.3, 0.3)
alpha1.true <- 0.12
alpha2.true <- 0.23
alpha3.true <- 0.34
kappa1.true <- 0.33
kappa2.true <- 0.11
kappa3.true <- 0.22
theta.true <- 0

# Make design matrix with single binary covariate
x_c <- rbinom(n, size = 1, prob = 0.7)
x_m <- cbind(1, x_c)

# Generate data
dat_ID <- simID(x1 = x_m, x2 = x_m, x3 = x_m,
                beta1.true = beta1.true, 
                beta2.true = beta2.true, 
                beta3.true = beta3.true,
                alpha1.true = alpha1.true, 
                alpha2.true = alpha2.true, 
                alpha3.true = alpha3.true,
                kappa1.true = kappa1.true, 
                kappa2.true = kappa2.true, 
                kappa3.true = kappa3.true,
                theta.true = theta.true,
                cens = c(240, 360))
dat_ID$x_c <- x_c
colnames(dat_ID)[1:4] <- c("R", "delta_R", "T", "delta_T")


# Plot --------------------------------------------------------------------

par(mfrow = c(1,2))
plot(survfit(Surv(R, delta_R) ~ x_c, data = dat_ID), mark.time = TRUE,
     col = hue_pal()(2), main = "Non-terminal", ylim = c(0, 1))
plot(survfit(Surv(T, delta_T) ~ x_c, data = dat_ID), mark.time = TRUE,
     col = hue_pal()(2), main = "Terminal", ylim = c(0, 1))
par(mfrow = c(1,1))

?SemiCompRisks::BayesSurv_HReg

# Data saving -------------------------------------------------------------

saveRDS(dat_ID, file = "./semicompstan/dat_ID.Rdata")



