#### 04. "What is this stuff called probability?
hdi_of_icdf <- function(name, width = .95, tol = 1e-8, ... ) {
  
  # Arguments:
  #   `name` is R's name for the inverse cumulative density function
  #   of the distribution.
  #   `width` is the desired mass of the HDI region.
  #   `tol` is passed to R's optimize function.
  # Return value:
  #   Highest density interval (HDI) limits in a vector.
  # Example of use: For determining HDI of a beta(30, 12) distribution, type
  #   `hdi_of_icdf(qbeta, shape1 = 30, shape2 = 12)`
  #   Notice that the parameters of the `name` must be explicitly stated;
  #   e.g., `hdi_of_icdf(qbeta, 30, 12)` does not work.
  # Adapted and corrected from Greg Snow's TeachingDemos package.
  
  incredible_mass <-  1.0 - width
  interval_width <- function(low_tail_prob, name, width, ...) {
    name(width + low_tail_prob, ...) - name(low_tail_prob, ...)
  }
  opt_info <- optimize(interval_width, c(0, incredible_mass), 
                       name = name, width = width, 
                       tol = tol, ...)
  hdi_lower_tail_prob <- opt_info$minimum
  
  return(c(name(hdi_lower_tail_prob, ...),
           name(width + hdi_lower_tail_prob, ...)))
  
}




library(tidyverse)
select <- dplyr::select

n       <- 500  # specify the total number of flips
p_heads <- 0.5  # specify underlying probability of heads

# Kruschke reported this was the seed he used at the top of page 94
set.seed(47405)

# here we use that seed to flip a coin n times and compute the running proportion of heads at each flip. 
# we generate a random sample of n flips (heads = 1, tails = 0)
d <-
  tibble(flip_sequence = sample(x = c(0, 1), 
                                prob = c(1 - p_heads, p_heads), 
                                size = n, 
                                replace = T)) %>% 
  mutate(n = 1:n,
         r = cumsum(flip_sequence)) %>% 
  mutate(run_prop = r / n)

end_prop <-
  d %>% 
  select(run_prop) %>% 
  slice(n()) %>% 
  round(digits = 3) %>% 
  pull()

d %>%
  filter(n < 1000) %>%  # this step cuts down on the time it takes to make the plot
  ggplot(aes(x = n, y = run_prop)) +
  geom_hline(yintercept = .5, color = "white") +
  geom_line(color = "grey50") +
  geom_point(color = "grey50", alpha = 1/4) +
  scale_x_log10(breaks = c(1, 2, 5, 10, 20, 50, 200, 500)) +
  coord_cartesian(xlim = c(1, 500),
                  ylim = c(0, 1)) +
  labs(title = "Running proportion of heads",
       subtitle = paste("Our end proportion =", end_prop),
       x = "Flip number",
       y = "Proportion of heads") +
  theme(panel.grid = element_blank())


HtWtDataGenerator <- function(n_subj, rndsd = NULL, male_prob = 0.50) {
  
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
  data_matrix <- matrix(0, nrow = n_subj, ncol = 3)
  colnames(data_matrix) <- c("male", "height", "weight")
  maleval <- 1; femaleval <- 0 # arbitrary coding values
  for (i in 1:n_subj) {
    # Flip coin to decide sex
    sex = sample(c(maleval, femaleval), size = 1, replace = TRUE,
                 prob = c(male_prob, 1 - male_prob))
    if (sex == maleval) {datum <- MASS::mvrnorm(n = 1, mu = Mmean, Sigma = Msigma)}
    if (sex == femaleval) {
      Fclust = sample(c(1, 2), size = 1, replace = TRUE, prob = c(prop1, prop2))
      if (Fclust == 1) {datum <- MASS::mvrnorm(n = 1, mu = Fmean1, Sigma = Fsigma1)}
      if (Fclust == 2) {datum <- MASS::mvrnorm(n = 1, mu = Fmean2, Sigma = Fsigma2)}
    }
    data_matrix[i, ] = c(sex, round(c(datum[1], exp(datum[2])), 1))
  }
  
  return(data_matrix)
  
} # end function


set.seed(4)
d <-
  HtWtDataGenerator(n_subj = 10000, male_prob = .5) %>%
  data.frame() %>%
  mutate(person = 1:n())

d %>%
  head()

d_bin <-
  d %>%
  mutate(bin = case_when(
    height < 51 ~ 51,
    between(height, 51, 53) ~ 53,
    between(height, 53, 55) ~ 55,
    between(height, 55, 57) ~ 57,
    between(height, 57, 59) ~ 59,
    between(height, 59, 61) ~ 61,
    between(height, 61, 63) ~ 63,
    between(height, 63, 65) ~ 65,
    between(height, 65, 67) ~ 67,
    between(height, 67, 69) ~ 69,
    between(height, 69, 71) ~ 71,
    between(height, 71, 73) ~ 73,
    between(height, 73, 75) ~ 75,
    between(height, 75, 77) ~ 77,
    between(height, 77, 79) ~ 79,
    between(height, 79, 81) ~ 71,
    between(height, 81, 83) ~ 83,
    height > 83 ~ 85)
  ) %>%
  group_by(bin) %>%
  summarise(n = n()) %>%
  mutate(height = bin - 1)

d %>%
  ggplot(aes(x = height, y = person)) +
  geom_point(size = 3/4, color = "grey67", alpha = 1/2) +
  geom_vline(xintercept = seq(from = 51, to = 83, by = 2),
             linetype = 3, color = "grey33") +
  geom_text(data = d_bin, 
            aes(y = 5000, label = n),
            size = 3.25) +
  scale_y_continuous(breaks = c(0, 5000, 10000)) +
  labs(title = "Total N = 10,000",
       x = "Height (inches)",
       y = "Person #") +
  theme(panel.grid = element_blank())



d %>%
  ggplot(aes(x = height)) +
  geom_histogram(aes(y = stat(density)),
                 binwidth = 2, fill = "grey67",
                 color = "grey92", size = 1/8) +
  scale_y_continuous("Probability density", breaks = c(0, .04, .08)) +
  xlab("Height (inches)") +
  coord_cartesian(xlim = c(51, 83)) +
  theme(panel.grid = element_blank()) 


d_bin <-
  d %>%
  mutate(bin = round(height, digits = 0)) %>%   
  group_by(bin) %>%
  summarise(n = n()) %>%
  mutate(height = bin - 0.5)

d %>%
  ggplot(aes(x = height, y = person)) +
  geom_point(size = 3/4, color = "grey67", alpha = 1/2) +
  geom_vline(xintercept = seq(from = 51, to = 83, by = 1),
             linetype = 3, color = "grey33") +
  geom_text(data = d_bin, 
            aes(y = 5000, label = n, angle = 90),
            size = 3.25) +
  scale_y_continuous(breaks = c(0, 5000, 10000)) +
  labs(title = "Total N = 10,000",
       x = "Height (inches)",
       y = "Person #") +
  theme(panel.grid = element_blank())

d %>%
  ggplot(aes(x = height)) +
  geom_histogram(aes(y = stat(density)), boundary = 0,
                 binwidth = 1, fill = "grey67",
                 color = "grey92", size = 1/8) +
  scale_y_continuous("Probability density", breaks = c(0, .04, .08)) +
  xlab("Height (inches)") +
  coord_cartesian(xlim = c(51, 83)) +
  theme(panel.grid = element_blank())


set.seed(4)
d <-
  tibble(height = rnorm(1e4, mean = 84, sd = .1)) %>%
  mutate(door = 1:n())

d %>%
  head()

library(santoku)

d_bin <-
  d %>% 
  mutate(bin = chop(height, 
                    breaks = seq(from = 83.6, to = 84.4, length.out = 31),
                    labels = seq(from = 83.6, to = 84.4, length.out = 31)[-1]))

head(d_bin)

d_bin <-
  d_bin %>% 
  mutate(bin = as.character(bin) %>% as.double()) %>% 
  group_by(bin) %>%
  summarise(n = n()) %>% 
  mutate(height = bin - (83.62667 - 83.6) / 2)

head(d_bin)

d %>%
  ggplot(aes(x = height, y = door)) +
  geom_point(size = 3/4, color = "grey67", alpha = 1/2) +
  geom_vline(xintercept = seq(from = 83.6, to = 84.4, length.out = 31),
             linetype = 3, color = "grey33") +
  geom_text(data = d_bin,
            aes(y = 5000, label = n, angle = 90),
            size = 3.25) +
  scale_y_continuous(breaks = c(0, 5000, 10000)) +
  labs(title = "Total N = 10,000",
       x = "Height (inches)",
       y = "Door #") +
  theme(panel.grid = element_blank())

d %>%
  ggplot(aes(x = height)) +
  geom_histogram(aes(y = stat(density)), boundary = 0,
                 binwidth = .025, fill = "grey67",
                 color = "grey92", size = 1/8) +
  scale_y_continuous("Probability density", breaks = 0:4) +
  xlab("Height (inches)") +
  coord_cartesian(xlim = c(83.6, 84.4)) +
  theme(panel.grid = element_blank())


tibble(x = seq(from = -.8, to = .8, by = .02)) %>% 
  mutate(p = dnorm(x, mean = 0, sd = .2)) %>% 
  ggplot(aes(x = x)) +
  geom_line(aes(y = p),
            color = "grey50", size = 1.25) +
  geom_linerange(aes(ymin = 0, ymax = p),
                 size = 1/3) +
  labs(title = "Normal probability density",
       subtitle = expression(paste(mu, " = 0 and ", sigma, " = 0.2")),
       y = "p(x)") +
  coord_cartesian(xlim = c(-.61, .61)) +
  theme(panel.grid = element_blank())

h <-
  hdi_of_icdf(name = qnorm,
              mean = 0,
              sd   = 1)

h


tibble(x = seq(from = -3.5, to = 3.5, by = .05)) %>% 
  mutate(d = dnorm(x, mean = 0, sd = 1)) %>% 
  
  ggplot(aes(x = x, y = d)) +
  geom_area(fill = "grey75") +
  geom_area(data = . %>% filter(x >= h[1] & x <= h[2]),
            fill = "grey50") +
  geom_line(data = tibble(x = c(h[1] + .02, h[2] - .02),
                          d = c(.059, .059)),
            arrow = arrow(length = unit(.2, "cm"), 
                          ends = "both", 
                          type = "closed"),
            color = "grey92") +
  annotate(geom = "text", x = 0, y = .09, 
           label = "95% HDI", color = "grey92") +
  xlim(-3.1, 3.1) +
  ylab("p(x)") +
  theme(panel.grid = element_blank())

h <-
  hdi_of_icdf(name = qbeta,
              shape1 = 15, 
              shape2 = 4)

h


tibble(x = seq(from = 0, to = 1, by = .01)) %>% 
  mutate(d = dbeta(x, shape1 = 15, shape2 = 4)) %>% 
  
  ggplot(aes(x = x, y = d)) +
  geom_area(fill = "grey75") +
  geom_area(data = . %>% filter(x >= h[1] & x <= h[2]),
            fill = "grey50") +
  geom_line(data = tibble(x = c(h[1] + .01, h[2] - .002),
                          d = c(.75, .75)),
            arrow = arrow(length = unit(.2, "cm"),
                          ends = "both",
                          type = "closed"),
            color = "grey92") +
  annotate(geom = "text", x = .8, y = 1.1, 
           label = "95% HDI", color = "grey92") +
  xlim(.4, 1) +
  ylab("p(x)") +
  theme(panel.grid = element_blank())


set.seed(4)
d <-
  tibble(x = c(rnorm(6e5, mean = 1.50, sd = .5),
               rnorm(4e5, mean = 4.75, sd = .5)))

glimpse(d)

library(tidybayes)

h <- 
  d %>% 
  mode_hdi()

h

dens <-
  d$x %>%
  density() %>%
  with(tibble(x, y))

head(dens)

ggplot(data = dens,
       aes(x = x, y = y)) +
  geom_area(fill = "grey75") +
  # note the use of `pull()`, which extracts the values, rather than return a tibble  
  geom_area(data = dens %>% filter(x > h[1, 2] %>% pull() & 
                                     x < h[1, 3] %>% pull()),
            fill = "grey50") +
  geom_area(data = dens %>% filter(x > h[2, 2] %>% pull() & 
                                     x < h[2, 3] %>% pull()),
            fill = "grey50") +
  geom_line(data = tibble(x = c(h[1, 2] %>% pull(), h[1, 3] %>% pull()),
                          y = c(.06, .06)),
            arrow = arrow(length = unit(.2,"cm"),
                          ends = "both",
                          type = "closed"),
            color = "grey92") +
  geom_line(data = tibble(x = c(h[2, 2] %>% pull(), h[2, 3] %>% pull()),
                          y = c(.06, .06)),
            arrow = arrow(length = unit(.2,"cm"),
                          ends = "both",
                          type = "closed"),
            color = "grey92") +
  annotate(geom = "text", x = c(1.5, 4.75), y = .1, 
           label = "95% HDI", color = "grey92") +
  scale_x_continuous(breaks = 0:6, limits = c(0, 6.3)) +
  scale_y_continuous("p(x)", breaks = c(0, .1, .2, .3, .4, .5)) +
  theme(panel.grid = element_blank())
