#### 06. Inferring a binomial probability via exact mathematical analysis
source('./scripts/setup.R')

theta <- .5
a     <- 3
b     <- 3

dbeta(theta, a, b)
beta(a, b)


library(tidyverse)

length <- 1e4

d <-
  crossing(shape1 = c(.1, 1:4),
           shape2 = c(.1, 1:4)) %>%
  expand(nesting(shape1, shape2),
         x = seq(from = 0, to = 1, length.out = length)) %>% 
  mutate(a     = str_c("a = ", shape1),
         b     = str_c("b = ", shape2),
         group = rep(1:length, each = 25))

head(d)

d %>% 
  ggplot(aes(x = x, group = group)) +
  
  geom_line(aes(y = dbeta(x, shape1 = shape1, shape2 = shape2)),
            color = "grey50", size = 1.25) +
  scale_x_continuous(expression(theta), breaks = c(0, .5, 1)) +
  coord_cartesian(ylim = c(0, 3)) +
  labs(title = "Examples of the beta distribution",
       y = expression(p(theta*"|"*a*", "*b))) +
  theme(panel.grid = element_blank()) +
  facet_grid(b ~ a)

d <-
  tibble(shape1 = c(5.6, 17.6, 5, 17),
         shape2 = c(1.4, 4.4, 2, 5)) %>% 
  mutate(a        = str_c("a = ", shape1),
         b        = str_c("b = ", shape2),
         kappa    = rep(c("kappa==7", "kappa==22"), times = 2),
         mu_omega = rep(c("mu==0.8", "omega==0.8"), each = 2)) %>% 
  mutate(kappa = factor(kappa, levels = c("kappa==7", "kappa==22")),
         label = str_c(a, ", ", b)) %>% 
  expand(nesting(shape1, shape2, a, b, label, kappa, mu_omega), 
         x = seq(from = 0, to = 1, length.out = length))

head(d)

d %>%
  ggplot(aes(x = x)) +
  geom_vline(xintercept = .8, color = "white") +
  geom_line(aes(y = dbeta(x, shape1 = shape1, shape2 = shape2)),
            color = "grey50", size = 1.25) +
  
  geom_text(data = . %>% group_by(label) %>% slice(1),
            aes(x = .025, y = 4.75, label = label),
            hjust = 0, size = 3) +
  
  scale_x_continuous(expression(theta), breaks = c(0, .8, 1)) +
  ylab(expression(p(theta*"|"*a*", "*b))) +
  coord_cartesian(ylim = c(0, 5)) +
  theme(panel.grid = element_blank()) +
  facet_grid(mu_omega ~ kappa, labeller = label_parsed)


betaAB(x = .25, y = 4, method = 'mean')

betaAB(x = .25, y = 4, method = 'mode')

betaAB(x = 0.5, y = .1, method = 'sd')


beta_param <- betaAB(x = .25, y = 4, method = 'mode')

beta_param$a
beta_param$b

hdi_of_icdf(name = qbeta,
            shape1 = 5,
            shape2 = 5,
            width  = .95)

hdi_of_icdf(name = qbeta,
            shape1 = 6,
            shape2 = 14)

n <- 10
z <- 1
a <- 5
b <- 5

(proportion_heads <- z / n)

(prior_mean <- a / (a + b))

(posterior_mean <- (z + a) / (n + a + b))


trial_data <- c(rep(0, 9), 1)

d <-
  tibble(theta = seq(from = 0, to = 1, length.out = 100)) %>% 
  mutate(`Prior (beta)`           = dbeta(theta, 
                                          shape1 = a, 
                                          shape2 = b),
         `Likelihood (Bernoulli)` = bernoulli_likelihood(theta = theta, 
                                                         data  = trial_data),
         `Posterior (beta)`       = dbeta(theta, 
                                          shape1 = 6, 
                                          shape2 = 14))

glimpse(d)

# the data for the in-plot lines
line <-
  tibble(theta      = c(.212 + .008, .788 - .008, .114 + .004, .497 - .005),
         value      = rep(c(.51, .66), each = 2),
         xintercept = c(.212, .788, .114, .497),
         key        = rep(c("Prior (beta)", "Posterior (beta)"), each = 2)) %>% 
  mutate(key = factor(key, levels = c("Prior (beta)", "Likelihood (Bernoulli)", "Posterior (beta)")))

# the data for the annotation
text <-
  tibble(theta = c(.5, .3),
         value = c(.8, 1.125),
         label = "95% HDI",
         key   = c("Prior (beta)", "Posterior (beta)")) %>% 
  mutate(key = factor(key, levels = c("Prior (beta)", "Likelihood (Bernoulli)", "Posterior (beta)")))


library(cowplot)

d %>% 
  gather(key, value, -theta) %>% 
  mutate(key = factor(key, levels = c("Prior (beta)", "Likelihood (Bernoulli)", "Posterior (beta)"))) %>% 
  
  ggplot(aes(x = theta, y = value, )) +
  # densities
  geom_area(fill = utastoolkit::utas_pal(6, 'new')[3]) +
  # dashed vertical lines
  geom_vline(data = line,
             aes(xintercept = xintercept), 
             linetype = 2, color = "white") +
  # arrows
  geom_line(data = line,
            arrow = arrow(length = unit(.15,"cm"), 
                          ends = "both", 
                          type = "closed"),
            color = "white") +
  # text
  geom_text(data = text,
            aes(label = label),
            color = "white") +
  labs(x = expression(theta),
       y = NULL) +
  facet_wrap(~ key, scales = "free_y", ncol = 1) +
  theme_cowplot()

beta_param <- 
  betaAB(x = 0.5, y = 200,
         method = 'mode')

# compute the corresponding HDIs
prior_hdi <-
  hdi_of_icdf(name = qbeta,
              shape1 = beta_param$a,
              shape2 = beta_param$b,
              width = .95)

# define the data
n <- 20
z <- 17

trial_data <- c(rep(0, times = n - z), rep(1, times = z))

# compute the HDIs for the posterior
post_hdi <-
  hdi_of_icdf(name = qbeta,
              shape1 = z + beta_param$a,
              shape2 = n - z + beta_param$b,
              width = .95)

# use the above to compute the prior, the likelihood, and the posterior
# densities using the grid approximation approach
d <-
  tibble(theta = seq(from = 0, to = 1, length.out = 1e3)) %>% 
  mutate(prior      = dbeta(theta, 
                            shape1 = beta_param$a, 
                            shape2 = beta_param$b),
         likelihood = bernoulli_likelihood(theta = theta, 
                                           data  = trial_data),
         posterior  = dbeta(theta, 
                            shape1 = z + beta_param$a, 
                            shape2 = n - z + beta_param$b))

# what have we done?
glimpse(d)

## Figure 6.4, left column
# prior
d %>% 
  ggplot(aes(x = theta, y = prior)) +
  geom_area(fill = utastoolkit::utas_pal(3)[3], alpha = 1/2) +
  geom_area(data = . %>% filter(theta > prior_hdi[1] & theta < prior_hdi[2]),
            fill = utastoolkit::utas_pal(3)[3]) +
  geom_segment(x = prior_hdi[1] + .005, xend = prior_hdi[2] - .005,
               y = 1.8, yend = 1.8,
               arrow = arrow(length = unit(.15,"cm"), 
                             ends = "both", 
                             type = "closed"),
               color = "white") +
  annotate(geom = "text", x = .5, y = 3.5, 
           label = "95% HDI") +
  labs(title = "Prior (beta)",
       x = expression(theta),
       y = expression(dbeta(theta*"|"*100*", "*100))) +
  coord_cartesian(ylim = c(0, 12)) +
  theme_cowplot()


d %>%   
  ggplot(aes(x = theta, y = likelihood)) +
  geom_area(fill = "steelblue") +
  labs(title = "Likelihood (Bernoulli)",
       x = expression(theta),
       y = expression(p(D*"|"*theta))) +
  theme_cowplot()


d %>% 
  ggplot(aes(x = theta, y = posterior)) +
  geom_area(fill = "steelblue", alpha = 1/2) +
  geom_area(data = . %>% filter(theta > post_hdi[1] & theta < post_hdi[2]),
            fill = "steelblue") +
  geom_segment(x = post_hdi[1] + .005, xend = post_hdi[2] - .005,
               y = 2, yend = 2,
               arrow = arrow(length = unit(.15, "cm"), 
                             ends = "both", 
                             type = "closed"),
               color = "white") +
  annotate(geom = "text", x = .532, y = 3.5, 
           label = "95% HDI") +
  labs(title = "Posterior (beta)",
       x = expression(theta),
       y = expression(dbeta(theta*"|"*117*", "*103))) +
  coord_cartesian(ylim = c(0, 12)) +
  theme_cowplot()


prior_hdi
post_hdi

(beta_param$a - 1) / (beta_param$a + beta_param$b - 2)

(z + beta_param$a - 1) / (z + beta_param$a + n - z + beta_param$b - 2)

# update the beta parameters for the prior
beta_param$a <- 18.25
beta_param$b <- 6.75

# update the HDIs
prior_hdi <-
  hdi_of_icdf(name = qbeta,
              shape1 = beta_param$a,
              shape2 = beta_param$b,
              width = .95)

post_hdi <-
  hdi_of_icdf(name = qbeta,
              shape1 = z + beta_param$a,
              shape2 = n - z + beta_param$b,
              width = .95)

# update the data
d <-
  d %>% 
  mutate(prior     = dbeta(theta, 
                           shape1 = beta_param$a, 
                           shape2 = beta_param$b),
         posterior = dbeta(theta, 
                           shape1 = z + beta_param$a, 
                           shape2 = n - z + beta_param$b))


## plot Figure 6.4, middle column!
# prior
d %>% 
  ggplot(aes(x = theta, y = prior)) +
  geom_area(fill = "steelblue", alpha = 1/2) +
  geom_area(data = . %>% filter(theta > prior_hdi[1] & theta < prior_hdi[2]),
            fill = "steelblue") +
  geom_segment(x = prior_hdi[1] + .005, xend = prior_hdi[2] - .005,
               y = 0.75, yend = 0.75,
               arrow = arrow(length = unit(.15,"cm"), 
                             ends = "both", 
                             type = "closed"),
               color = "white") +
  annotate(geom = "text", x = .75, y = 1.5, 
           label = "95% HDI", color = "white") +
  labs(title = "Prior (beta)",
       x = expression(theta),
       y = expression(dbeta(theta*"|"*18.25*", "*6.75))) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_cowplot()


# likelihood, which is the same as the last time
d %>%   
  ggplot(aes(x = theta, y = likelihood)) +
  geom_area(fill = "steelblue") +
  labs(title = "Likelihood (Bernoulli)",
       x = expression(theta),
       y = expression(p(D*"|"*theta))) +
  theme_cowplot()

# posterior
d %>% 
  ggplot(aes(x = theta, y = posterior)) +
  geom_area(fill = "steelblue", alpha = 1/2) +
  geom_area(data = . %>% filter(theta > post_hdi[1] & theta < post_hdi[2]),
            fill = "steelblue") +
  geom_segment(x = post_hdi[1] + .005, xend = post_hdi[2] - .005,
               y = 1, yend = 1,
               arrow = arrow(length = unit(.15, "cm"), 
                             ends = "both", 
                             type = "closed"),
               color = "white") +
  annotate(geom = "text", x = .797, y = 2, 
           label = "95% HDI", color = "white") +
  labs(title = "Posterior (beta)",
       x = expression(theta),
       y = expression(dbeta(theta*"|"*35.25*", "*9.75))) +
  coord_cartesian(ylim = c(0, 7)) +
  theme_cowplot()


prior_hdi
post_hdi
(beta_param$a - 1) / (beta_param$a + beta_param$b - 2)
(z + beta_param$a - 1) / (z + beta_param$a + n - z + beta_param$b - 2)
# update beta_param
beta_param$a <- 1
beta_param$b <- 1

# update the HDIs
prior_hdi <-
  hdi_of_icdf(name = qbeta,
              shape1 = beta_param$a,
              shape2 = beta_param$b,
              width = .95)

post_hdi <-
  hdi_of_icdf(name = qbeta,
              shape1 = z + beta_param$a,
              shape2 = n - z + beta_param$b,
              width = .95)

# update the data
d <-
  d %>% 
  mutate(prior     = dbeta(theta, 
                           shape1 = beta_param$a, 
                           shape2 = beta_param$b),
         posterior = dbeta(theta, 
                           shape1 = z + beta_param$a, 
                           shape2 = n - z + beta_param$b))


## plot Figure 6.4, rightmost column!
# prior
d %>% 
  ggplot(aes(x = theta, y = prior)) +
  geom_area(fill = "steelblue") +
  labs(title = "Prior (beta)",
       x = expression(theta),
       y = expression(dbeta(theta*"|"*1*", "*1))) +
  coord_cartesian(ylim = c(0, 5)) +
  theme_cowplot()

# likelihood, which is the same as the last two examples
d %>%   
  ggplot(aes(x = theta, y = likelihood)) +
  geom_area(fill = "steelblue") +
  labs(title = "Likelihood (Bernoulli)",
       x = expression(theta),
       y = expression(p(D*"|"*theta))) +
  theme_cowplot()

# posterior
d %>% 
  ggplot(aes(x = theta, y = posterior)) +
  geom_area(fill = "steelblue", alpha = 1/2) +
  geom_area(data = . %>% filter(theta > post_hdi[1] & theta < post_hdi[2]),
            fill = "steelblue") +
  geom_segment(x = post_hdi[1] + .005, xend = post_hdi[2] - .005,
               y = 0.8, yend = 0.8,
               arrow = arrow(length = unit(.15, "cm"), 
                             ends = "both", 
                             type = "closed"),
               color = "white") +
  annotate(geom = "text", x = (post_hdi[1] + post_hdi[2]) / 2, y = 1.5, 
           label = "95% HDI", color = "white") +
  labs(title = "Posterior (beta)",
       x = expression(theta),
       y = expression(dbeta(theta*"|"*18*", "*4))) +
  coord_cartesian(ylim = c(0, 5)) +
  theme_cowplot()


post_hdi
(z + beta_param$a - 1) / (z + beta_param$a + n - z + beta_param$b - 2)

# Fine teeth for Theta
theta <- seq(0, 1, length = 1000)

# Two triangular peaks on a small non-zero floor
p_theta <-
  c(rep(1, 200), 
    seq(1, 100, length = 50), 
    seq(100, 1, length = 50), 
    rep(1, 200)) %>% 
  rep(., times = 2)

# Make p_theta sum to 1.0
p_theta <- p_theta / sum(p_theta)

Data <- c(rep(0, 13), rep(1, 14))

BernGrid(theta, p_theta, Data, plotType = "Bars",
         showCentTend = "None", showHDI = FALSE, showpD = FALSE)


# we need these to compute the likelihood
n <- 27
z <- 14

trial_data <- c(rep(0, times = n - z), rep(1, times = z))        # (i.e., Data)

d <-
  tibble(theta = seq(from = 0, to = 1, length.out = 1000),       # (i.e., Theta)
         Prior = c(rep(1, 200),                                  # (i.e., pTheta)
                   seq(1, 100, length = 50), 
                   seq(100, 1, length = 50), 
                   rep(1, 200)) %>% 
           rep(., times = 2)) %>% 
  mutate(Prior      = Prior / sum(Prior),
         Likelihood = bernoulli_likelihood(theta = theta,        # (i.e., pDataGivenTheta)
                                           data  = trial_data)) %>%
  mutate(evidence = sum(Likelihood * Prior)) %>%                 # (i.e., pData)
  mutate(Posterior = Likelihood * Prior / evidence)              # (i.e., pThetaGivenData)

glimpse(d)

# prior
(p1 <-
    d %>% 
    ggplot(aes(x = theta, y = Prior)) +
    geom_area(fill = "steelblue") +
    labs(title = "Prior",
         x = expression(theta),
         y = expression(p(theta))) +
    theme_cowplot()
)

# likelihood
(p2 <-
    d %>% 
    ggplot(aes(x = theta, y = Likelihood)) +
    geom_area(fill = "steelblue") +
    labs(title = "Likelihood",
         x = expression(theta),
         y = expression(p(D*"|"*theta))) +
    theme_cowplot()
)
# posterior
(p3 <-
    d %>% 
    ggplot(aes(x = theta, y = Posterior)) +
    geom_area(fill = "steelblue") +
    labs(title = "Posterior",
         x = expression(theta),
         y = expression(p(theta*"|"*D))) +
    theme_cowplot()
)

library(patchwork)

p1 / p2 / p3
