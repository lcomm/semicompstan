library("SemiCompRisks")
library("mvnfast")
library("scales")
library("plyr")
library("dplyr")
library("ggplot2")
library("ggsci")
library("gridExtra")
library("grid")
library("ggpubr")

set.seed(42)


# Simulate data -----------------------------------------------------------

# Function to simulate data
simulate_RI <- function(n, 
                        theta.true, 
                        alpha1.true, alpha2.true, alpha3.true,
                        beta1.true, beta2.true, beta3.true,
                        kappa1.true, kappa2.true, kappa3.true,
                        cens) {
  
  # Cens -> vector
  if (length(cens) != 2) {
    cens <- c(cens, cens + 1)
  }
  
  # Frailty parameters
  log_gamma_i <- log(rgamma(n, shape = 1/theta.true, rate = 1/theta.true))
  
  # Simulate potential confounder
  x_c <- rnorm(n, mean = 0, sd = 1)
  a <- rep(0:1, each = n)
  x1 <- x2 <- x3 <- cbind(c(x_c, x_c), a, c(log_gamma_i, log_gamma_i))
  
  # Simulate, adding one to each set of parameters
  sim_dat <- simID(cluster = NULL, x1, x2, x3, 
                   alpha1.true = alpha1.true,
                   alpha2.true = alpha2.true,
                   alpha3.true = alpha3.true,
                   beta1.true = c(beta1.true, 1),
                   beta2.true = c(beta2.true, 1),
                   beta3.true = c(beta3.true, 1),
                   kappa1.true = kappa1.true,
                   kappa2.true = kappa2.true,
                   kappa3.true = kappa3.true,
                   theta.true, SigmaV.true = NULL, cens)
  
  # Rename and add treatment column
  colnames(sim_dat) <- c("R", "delta_R", "T", "delta_T")
  sim_dat$ID <- c(1:n, 1:n)
  sim_dat$z <- factor(rep(0:1, each = n))
  
  # Return created data
  return(sim_dat)
}

# True parameters
alpha1.true <- 0.2
alpha2.true <- 0.05
alpha3.true <- 0.15

beta3.true <- beta2.true <- beta1.true <- c(0, -0.4)
kappa3.true <- kappa2.true <- kappa1.true <- 0.1
theta.true <- 0.5
maxt <- 240
cens <- 10^6
n <- 10^3

dat <- simulate_RI(n, theta.true, 
                   alpha1.true, alpha2.true, alpha3.true,
                   beta1.true, beta2.true, beta3.true,
                   kappa1.true, kappa2.true, kappa3.true,
                   cens)


# Plotting functions ------------------------------------------------------

grid_arrange_shared_legend <- function(..., nrow = 1, ncol = length(list(...)), position = c("bottom", "right")) {
  
  plots <- list(...)
  position <- match.arg(position)
  g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
  legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
  lheight <- sum(legend$height)
  lwidth <- sum(legend$width)
  gl <- lapply(plots, function(x) x + theme(legend.position = "none"))
  gl <- c(gl, nrow = nrow, ncol = ncol)
  
  combined <- switch(position,
                     "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                            legend,
                                            ncol = 1,
                                            heights = unit.c(unit(1, "npc") - lheight, lheight)),
                     "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                           legend,
                                           ncol = 2,
                                           widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
  grid.newpage()
  grid.draw(combined)
  
}


# Process simulated data --------------------------------------------------


# Get principal state
get_state <- function(T0, T1, at_time) {
  if ((T0 > at_time) & (T1 > at_time)) {
    state_t <- "AA"
  } else if ((T0 > at_time) & (T1 < at_time)) {
    state_t <- "TK"
  } else if ((T0 < at_time) & (T1 > at_time)) {
    state_t <- "CK"
  } else if ((T0 < at_time) & (T1 < at_time)) {
    state_t <- "DD"
  } else {
    state_t <- "undef"
  }
  return(state_t)
}

dat$World <- factor(dat$z, levels = c("1", "0"), 
                    labels = c("Treated (z=1)", "Control (z=0)"))

rbarnum <- maxt + 30
maxtplot <- rbarnum + 5

pctRbar <- group_by(dat, World) %>% 
                       summarise(mean(delta_R == 0, na.rm = TRUE))
pctRbar$rbarnum <- rbarnum

dat$Rtoplot <- dat$R
dat$Rtoplot[dat$delta_R == 0] <- rbarnum * 2
rplot <- ggplot(dat, aes(Rtoplot, fill = World, colour = World)) +
    geom_density(alpha = 0.2, trim = maxt) + 
    ylab("Density") + 
    ggtitle("Rehospitalization")
showRbar <- TRUE
if (showRbar) {
  rplot <- rplot  +
    scale_x_continuous("Rehospitalization time",
                       limits = c(0, maxtplot),
                       breaks = c(seq(0, maxt, by = 30), rbarnum),
                       labels = c(seq(0, maxt, by = 30), expression(bold(bar(R))))) +
    annotate("segment", color = hue_pal()(2)[which.max(pctRbar[[2]])], 
             size = 1.5,
             x = rbarnum, xend = rbarnum, 
             y = 0, yend = max(pctRbar[[2]])) +
    annotate("segment", color = hue_pal()(2)[which.min(pctRbar[[2]])], 
             size = 1.5,
             x = rbarnum, xend = rbarnum, 
             y = 0, yend = min(pctRbar[[2]]))
} else {
  rplot <- rplot + 
    scale_x_continuous("Rehospitalization time",
                       limits = c(0, maxtplot),
                       breaks = seq(0, maxt, by = 30),
                       labels = seq(0, maxt, by = 30))
}
  
tplot <- ggplot(dat, aes(T, fill = World, colour = World)) +
  geom_density(alpha = 0.2) +
  ylab("Density") + 
  ggtitle("Death") +
  scale_x_continuous("Death time",
                     limits = c(0, maxt),
                     breaks = seq(0, maxt, by = 30),
                     labels = seq(0, maxt, by = 30))

suppressWarnings(grid_arrange_shared_legend(rplot, tplot, nrow = 1, ncol = 2))



# Composition stream plots ------------------------------------------------

nxts <- 10^3
xts <- seq(0, maxt, length.out = nxts)

states <- diffs <- matrix(NA, nrow = n, ncol = nxts)

for (id in 1:n) {
  T0 <- dat$T[dat$ID == id & dat$z == 0]
  T1 <- dat$T[dat$ID == id & dat$z == 1]
  R0 <- dat$R[dat$ID == id & dat$z == 0]
  R1 <- dat$R[dat$ID == id & dat$z == 1]
  
  states[id, ] <- sapply(xts, FUN = get_state, T0 = T0, T1 = T1)
  diffs[id, ]  <- (R1 < xts) - (R0 < xts)
  
}


# Data frame for composition streamgraph
comp_dat0 <- data.frame(time = xts)
pstates <- c("AA", "CK", "TK", "DD")
for (pstate in pstates) {
  comp_dat0[[pstate]] <- colMeans(states == pstate)
}

comp_dat <- data.frame(Time = rep(xts, times = length(pstates)),
                       State = rep(pstates, each = nxts))

comp_dat$Proportion <- c(comp_dat0[["AA"]], comp_dat0[["CK"]],
                         comp_dat0[["TK"]], comp_dat0[["DD"]])
comp_dat$State <- factor(comp_dat$State, levels = rev(pstates), 
                         ordered = TRUE)

cplot <- ggplot(comp_dat, aes(x = Time, y = Proportion, fill = State)) + 
    geom_area(alpha = 0.6, color = "black") + 
    scale_x_continuous(breaks = seq(0, maxt, by = 30)) +
    ggtitle("Principal state composition over time") + 
    scale_fill_discrete(breaks = rev(levels(comp_dat$State))) + 
    theme(legend.position = "bottom", legend.key.size = unit(0.6,"line"))


# Causal effect plot ------------------------------------------------------

diffs_AA <- diffs
diffs_AA[states != "AA"] <- NA

ce_dat <- data.frame(Time = xts,
                     alpha_tt = colMeans(diffs_AA, na.rm = TRUE))

aplot <- ggplot(ce_dat, aes(x = Time, y = alpha_tt)) +
    geom_smooth(method = "loess", se = FALSE) + 
    ylab(expression(alpha(r,t)~with~r==t)) + 
    scale_x_continuous(breaks = seq(0, maxt, by = 30)) +
    ggtitle("Causal effect on rehospitalization among Always Alive at t") +
    labs(subtitle = expression(alpha(r,t) == P(R[1]<r~'|'~T[0] > t, T[1] > t)-
                               P(R[0]<r~'|'~T[0] > t, T[1] > t)))
         
grid.arrange(cplot, aplot, nrow = 2, ncol = 1)

