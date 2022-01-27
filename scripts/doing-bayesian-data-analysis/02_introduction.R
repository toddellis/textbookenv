#### HIGHLIGHTS ----------------------------------------------------------------
### HDI in ggplot!!!! ----------------------------------------------------------
# stat_histinterval(point_interval = mode_hdi, 
#                   .width = .95,
#                   fill = "grey67",
#                   slab_color = "grey92",
#                   breaks = 40,
#                   slab_size = .2, 
#                   outline_bars = T)
## -----------------------------------------------------------------------------


library(tidyverse)

d <-
  crossing(iteration = 1:3,
           stage     = factor(c("Prior", "Posterior"),
                              levels = c("Prior", "Posterior"))) %>% 
  expand(nesting(iteration, stage),
         Possibilities = LETTERS[1:4]) %>% 
  mutate(Credibility   = c(rep(.25, times = 4),
                           0, rep(1/3, times = 3),
                           0, rep(1/3, times = 3),
                           rep(c(0, .5), each = 2),
                           rep(c(0, .5), each = 2),
                           rep(0, times = 3), 1))

text <-
  tibble(Possibilities = "B",
         Credibility   = .75,
         label         = str_c(LETTERS[1:3], " is\nimpossible"),
         iteration     = 1:3,
         stage         = factor("Posterior", 
                                levels = c("Prior", "Posterior")))
arrow <-
  tibble(Possibilities = LETTERS[1:3],
         iteration     = 1:3) %>% 
  expand(nesting(Possibilities, iteration),
         Credibility = c(0.6, 0.01)) %>% 
  mutate(stage = factor("Posterior", levels = c("Prior", "Posterior")))

d %>%
  ggplot(aes(x = Possibilities, y = Credibility)) +
  geom_col(color = "grey30", fill = "grey30") +
  # annotation in the bottom row
  geom_text(data = text, 
            aes(label = label)) +
  # arrows in the bottom row
  geom_line(data = arrow,
            arrow = arrow(length = unit(0.30, "cm"), 
                          ends = "first", type = "closed")) +
  facet_grid(stage ~ iteration) +
  theme(axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        strip.text.x = element_blank())

crossing(stage         = factor(c("Prior", "Posterior"),
                                levels = c("Prior", "Posterior")),
         Possibilities = LETTERS[1:4]) %>% 
  mutate(Credibility = c(rep(0.25, times = 4),
                         rep(0,    times = 3), 
                         1)) %>%
  
  # plot!
  ggplot(aes(x = Possibilities, y = Credibility)) +
  geom_col(color = "grey30", fill = "grey30") +
  # annotation in the bottom panel
  geom_text(data = tibble(
    Possibilities = "B",
    Credibility   = .8,
    label         = "D is\nresponsible",
    stage         = factor("Posterior", levels = c("Prior", "Posterior"))
  ), aes(label = label)
  ) +
  # the arrow
  geom_line(data = tibble(
    Possibilities = LETTERS[c(4, 4)],
    Credibility   = c(.25, .99),
    stage         = factor("Posterior", levels = c("Prior", "Posterior"))
  ),
  arrow = arrow(length = unit(0.30, "cm"), 
                ends = "last", type = "closed"),
  color = "grey92") +
  facet_wrap(~ stage, ncol = 1) +
  theme(axis.ticks.x = element_blank(),
        panel.grid = element_blank())


tibble(mu = 1:4,
       p  = .25) %>% 
  expand(nesting(mu, p), 
         x = seq(from = -2, to = 6, by = .1)) %>% 
  mutate(density = dnorm(x, mean = mu, sd = 1.2)) %>% 
  mutate(d_max = max(density)) %>% 
  mutate(rescale = p / d_max) %>% 
  mutate(density = density * rescale) %>% 
  # plot!
  ggplot(aes(x = x)) +
  geom_col(data = . %>% distinct(mu, p),
           aes(x = mu, y = p),
           fill = "grey67", width = 1/3) +
  geom_line(aes(y = density, group = mu)) +
  scale_x_continuous(breaks = 1:4) +
  scale_y_continuous(breaks = 0:5 / 5) +
  coord_cartesian(xlim = c(0, 5),
                  ylim = c(0, 1)) +
  labs(title = "Prior",
       x = "Possibilities", 
       y = "Credibility") +
  theme(axis.ticks.x = element_blank(),
        panel.grid = element_blank())


tibble(mu = 1:4,
       p  = c(.1, .58, .3, .02)) %>% 
  expand(nesting(mu, p), 
         x = seq(from = -2, to = 6, by = .1)) %>% 
  mutate(density = dnorm(x, mean = mu, sd = 1.2)) %>% 
  mutate(d_max = max(density)) %>% 
  mutate(rescale = p / d_max) %>% 
  mutate(density = density * rescale) %>% 
  
  # plot!
  ggplot() +
  geom_col(data = . %>% distinct(mu, p),
           aes(x = mu, y = p),
           fill = "grey67", width = 1/3) +
  geom_line(aes(x = x, y = density, group = mu)) +
  geom_point(data = tibble(x = c(1.75, 2.25, 2.75), y = 0),
             aes(x = x, y = y),
             size = 3, color = "grey33", alpha = 3/4) +
  scale_x_continuous(breaks = 1:4) +
  scale_y_continuous(breaks = 0:5 / 5) +
  coord_cartesian(xlim = c(0, 5),
                  ylim = c(0, 1)) +
  labs(title = "Posterior",
       x = "Possibilities", 
       y = "Credibility") +
  theme(axis.ticks.x = element_blank(),
        panel.grid = element_blank())

# set the seed to make the simulation reproducible
set.seed(2)
# simulate the data with `rnorm()`
d <- tibble(x = rnorm(2000, mean = 10, sd = 5))

# plot!
ggplot(data = d, aes(x = x)) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 1, fill = "grey67", 
                 color = "grey92", size = 1/10) +
  geom_line(data = tibble(x = seq(from = -6, to = 26, by = .01)),
            aes(x = x, y = dnorm(x, mean = 10, sd = 5)),
            color = "grey33") +
  coord_cartesian(xlim = c(-5, 25)) +
  labs(subtitle = "The candidate normal distribution\nhas a mean of 10 and SD of 5.",
       x = "Data Values", 
       y = "Data Probability") +
  theme(panel.grid = element_blank())

ggplot(data = d, aes(x = x)) +
  geom_histogram(aes(y = ..density..),
                 binwidth = 1, fill = "grey67",
                 color = "grey92", size = 1/8) +
  stat_function(fun = dnorm, n = 101, args = list(mean = 8, sd = 6),
                color = "grey33", linetype = 2) +
  coord_cartesian(xlim = c(-5, 25)) +
  labs(subtitle = "The candidate normal distribution\nhas a mean of 8 and SD of 6.",
       x = "Data Values", 
       y = "Data Probability") +
  theme(panel.grid = element_blank())


HtWtDataGenerator <- function(nSubj, rndsd = NULL, maleProb = 0.50) {
  # Random height, weight generator for males and females. Uses parameters from
  # Brainard, J. & Burmaster, D. E. (1992). Bivariate distributions for height and
  # weight of men and women in the United States. Risk Analysis, 12(2), 267-275.
  # Kruschke, J. K. (2011). Doing Bayesian data analysis:
  # A Tutorial with R and BUGS. Academic Press / Elsevier.
  # Kruschke, J. K. (2014). Doing Bayesian data analysis, 2nd Edition:
  # A Tutorial with R, JAGS and Stan. Academic Press / Elsevier.
  
  # require(MASS)
  
  # Specify parameters of multivariate normal (MVN) distributions.
  # Men:
  HtMmu   <- 69.18
  HtMsd   <- 2.87
  lnWtMmu <- 5.14
  lnWtMsd <- 0.17
  Mrho    <- 0.42
  Mmean   <- c(HtMmu, lnWtMmu)
  Msigma  <- matrix(c(HtMsd^2, Mrho * HtMsd * lnWtMsd,
                      Mrho * HtMsd * lnWtMsd, lnWtMsd^2), nrow = 2)
  # Women cluster 1:
  HtFmu1   <- 63.11
  HtFsd1   <- 2.76
  lnWtFmu1 <- 5.06
  lnWtFsd1 <- 0.24
  Frho1    <- 0.41
  prop1    <- 0.46
  Fmean1   <- c(HtFmu1, lnWtFmu1)
  Fsigma1  <- matrix(c(HtFsd1^2, Frho1 * HtFsd1 * lnWtFsd1,
                       Frho1 * HtFsd1 * lnWtFsd1, lnWtFsd1^2), nrow = 2)
  # Women cluster 2:
  HtFmu2   <- 64.36
  HtFsd2   <- 2.49
  lnWtFmu2 <- 4.86
  lnWtFsd2 <- 0.14
  Frho2    <- 0.44
  prop2    <- 1 - prop1
  Fmean2   <- c(HtFmu2, lnWtFmu2)
  Fsigma2  <- matrix(c(HtFsd2^2, Frho2 * HtFsd2 * lnWtFsd2,
                       Frho2 * HtFsd2 * lnWtFsd2, lnWtFsd2^2), nrow = 2)
  
  # Randomly generate data values from those MVN distributions.
  if (!is.null(rndsd)) {set.seed(rndsd)}
  datamatrix <- matrix(0, nrow = nSubj, ncol = 3)
  colnames(datamatrix) <- c("male", "height", "weight")
  maleval <- 1; femaleval <- 0 # arbitrary coding values
  for (i in 1:nSubj)  {
    # Flip coin to decide sex
    sex <- sample(c(maleval, femaleval), size = 1, replace = TRUE,
                  prob = c(maleProb, 1 - maleProb))
    if (sex == maleval) {datum = MASS::mvrnorm(n = 1, mu = Mmean, Sigma = Msigma)}
    if (sex == femaleval) {
      Fclust = sample(c(1, 2), size = 1, replace = TRUE, prob = c(prop1, prop2))
      if (Fclust == 1) {datum = MASS::mvrnorm(n = 1, mu = Fmean1, Sigma = Fsigma1)}
      if (Fclust == 2) {datum = MASS::mvrnorm(n = 1, mu = Fmean2, Sigma = Fsigma2)}
    }
    datamatrix[i, ] = c(sex, round(c(datum[1], exp(datum[2])), 1))
  }
  
  return(datamatrix)
} # end function

# set your seed to make the data generation reproducible
set.seed(2)

d <-
  HtWtDataGenerator(nSubj = 57, maleProb = .5) %>%
  as_tibble()

d %>%
  head()


library(brms)

fit2.1 <- 
  brm(data = d, 
      family = gaussian,
      weight ~ 1 + height,
      prior = c(prior(normal(0, 100), class = Intercept),
                prior(normal(0, 100), class = b),
                prior(cauchy(0, 10),  class = sigma)),
      chains = 4, cores = 4, iter = 2000, warmup = 1000,
      # file = "fits/fit02.01", ## This file doesn't exist without downloading it but doesn't seem necessary?? or does this *output* the file?
      seed = 2)
# extract the posterior draws
post <- posterior_samples(fit2.1)

# this will streamline some of the code, below
n_lines <- 150

# plot!
d %>% 
  ggplot(aes(x = height, y = weight)) +
  geom_abline(intercept = post[1:n_lines, 1], 
              slope     = post[1:n_lines, 2],
              color = "grey50", size = 1/4, alpha = .3) +
  geom_point(shape = 1) +
  # the `eval(substitute(paste()))` trick came from: https://www.r-bloggers.com/value-of-an-r-object-in-an-expression/
  labs(subtitle = eval(substitute(paste("Data with", n_lines, "credible regression lines"))),
       x = "Height in inches",
       y = "Weight in pounds") +
  coord_cartesian(xlim = c(55, 80),
                  ylim = c(50, 250)) +
  theme(panel.grid = element_blank())

library(tidybayes)

post %>% 
  ggplot(aes(x = b_height, y = 0)) +
  stat_histinterval(point_interval = mode_hdi, .width = .95,
                    fill = "grey67", slab_color = "grey92",
                    breaks = 40, slab_size = .2, outline_bars = T) +
  scale_y_continuous(NULL, breaks = NULL) +
  coord_cartesian(xlim = c(0, 8)) +
  labs(title = "The posterior distribution",
       subtitle = "The mode and 95% HPD intervals are\nthe dot and horizontal line at the bottom.",
       x = expression(beta[1]~(slope))) +
  theme(panel.grid = element_blank())

nd <- tibble(height = seq(from = 53, to = 81, length.out = 20))

predict(fit2.1, newdata = nd) %>%
  data.frame() %>%
  bind_cols(nd) %>% 
  
  ggplot(aes(x = height)) +
  geom_pointrange(aes(y = Estimate, ymin = Q2.5, ymax = Q97.5),
                  color = "grey67", shape = 20) +
  geom_point(data = d, 
             aes(y = weight),
             alpha = 2/3) +
  labs(subtitle = "Data with the percentile-based 95% intervals and\nthe means of the posterior predictions",
       x = "Height in inches",
       y = "Weight in inches") +
  theme(panel.grid = element_blank())


nd <- tibble(height = seq(from = 53, to = 81, length.out = 30))

predict(fit2.1, newdata = nd) %>%
  data.frame() %>%
  bind_cols(nd) %>% 
  
  ggplot(aes(x = height)) +
  geom_ribbon(aes(ymin = Q2.5, ymax = Q97.5), fill = "grey75") +
  geom_line(aes(y = Estimate),
            color = "grey92") +
  geom_point(data =  d, 
             aes(y = weight),
             alpha = 2/3) +
  labs(subtitle = "Data with the percentile-based 95% intervals and\nthe means of the posterior predictions",
       x = "Height in inches",
       y = "Weight in inches") +
  theme(panel.grid = element_blank())