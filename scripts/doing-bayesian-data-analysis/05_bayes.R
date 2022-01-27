#### 05. BAyes' Rule
bernoulli <- function(theta, y) {
  return(theta^y * (1 - theta)^(1 - y))
}

bernoulli_likelihood <- function(theta, data) {
  
  # `theta` = success probability parameter ranging from 0 to 1
  # `data` = the vector of data (i.e., a series of 0s and 1s)
  n   <- length(data)
  
  return(theta^sum(data) * (1 - theta)^(n - sum(data)))
  
}

library(tidyverse)

tibble(theta = seq(from = 0, to = 1, by = .1),
       prior = c(seq(from = 0, to = .2, length.out = 6),
                 seq(from = .16, to = 0, length.out = 5))) %>%
  
  ggplot(aes(x = theta, ymin = -0.0005, ymax = prior)) +
  geom_linerange(size = 4, color = "grey50") +
  scale_x_continuous(expression(theta), breaks = seq(from = 0, to = 1, by = .2)) +
  labs(title = "Prior",
       y = expression(p(theta))) +
  theme(panel.grid = element_blank())

tibble(prior = c(seq(from = 0,   to = .2, length.out = 6),
                 seq(from = .16, to = 0,  length.out = 5))) %>% 
  summarise(s = sum(prior))

theta_sequence <- seq(from = 0, to = 1, by = .1)

bernoulli(theta = theta_sequence, y = 1)

tibble(x = theta_sequence) %>%
  mutate(likelihood = bernoulli(theta = theta_sequence, y = 1)) %>% 
  
  ggplot(aes(x = x, y = likelihood)) +
  geom_col(width = .025, color = "grey50", fill = "grey50") +
  scale_x_continuous(expression(theta), breaks = seq(from = 0, to = 1, by = .2)) +
  labs(title = "Likelihood",
       y = expression(p(D*'|'*theta))) +
  theme(panel.grid = element_blank())


tibble(theta = theta_sequence,
       prior = c(seq(from = 0,   to = .2, length.out = 6),
                 seq(from = .16, to = 0,  length.out = 5))) %>%
  mutate(likelihood = bernoulli(theta = theta_sequence, y = 1)) %>% 
  mutate(marginal_likelihood = sum(prior * likelihood)) %>% 
  mutate(posterior = (prior * likelihood) / marginal_likelihood) %>%
  
  ggplot(aes(x = theta, y = posterior)) +
  geom_col(width = .025, color = "grey50", fill = "grey50") +
  scale_x_continuous(expression(theta), breaks = seq(from = 0, to = 1, by = .2)) +
  labs(title = "Posterior",
       y = expression(p(theta*'|'*D))) +
  theme(panel.grid = element_blank())

small_data <- rep(0:1, times = c(3, 1))

tibble(theta =   seq(from = 0,     to = 1, by = .001),
       Prior = c(seq(from = 0,     to = 1, length.out = 501),
                 seq(from = 0.998, to = 0, length.out = 500))) %>% 
  mutate(Prior      = Prior / sum(Prior),
         Likelihood = bernoulli_likelihood(theta = theta,
                                           data  = small_data)) %>% 
  mutate(marginal_likelihood = sum(Prior * Likelihood)) %>% 
  mutate(Posterior = (Prior * Likelihood) / marginal_likelihood) %>% 
  select(theta, Prior, Likelihood, Posterior) %>% 
  gather(key, value, -theta) %>% 
  mutate(key = factor(key, levels = c("Prior", "Likelihood", "Posterior"))) %>% 
  
  ggplot(aes(x = theta, y = value)) +
  geom_area(fill = "grey67") +
  scale_x_continuous(expression(theta), breaks = seq(from = 0, to = 1, by = .2)) +
  ylab("probability density") +
  theme(panel.grid = element_blank()) +
  facet_wrap(~ key, scales = "free_y", ncol = 1)

large_data <- rep(0:1, times = c(30, 10))

tibble(theta = seq(from = 0, to = 1, by = .001),
       Prior = c(seq(from = 0, to = 1, length.out = 501),
                 seq(from = 0.998, to = 0, length.out = 500))) %>% 
  mutate(Prior      = Prior / sum(Prior),
         Likelihood = bernoulli_likelihood(theta = theta, data = large_data)) %>% 
  mutate(marginal_likelihood = sum(Prior * Likelihood)) %>% 
  mutate(Posterior = (Prior * Likelihood) / marginal_likelihood) %>% 
  select(theta, Prior, Likelihood, Posterior) %>% 
  gather(key, value, -theta) %>% 
  mutate(key = factor(key, levels = c("Prior", "Likelihood", "Posterior"))) %>% 
  
  ggplot(aes(x = theta, y = value)) +
  geom_area(fill = "grey67") +
  scale_x_continuous(expression(theta), breaks = seq(from = 0, to = 1, by = .2)) +
  ylab("probability density") +
  theme(panel.grid = element_blank()) +
  facet_wrap(~ key, scales = "free_y", ncol = 1)


small_data <- rep(0:1, times = c(3, 1))

tibble(theta = seq(from = 0, to = 1, by = .001),
       Prior = c(seq(from = 0, to = 1, length.out = 501),
                 seq(from = 0.998, to = 0, length.out = 500))) %>% 
  # here's the important line of code
  mutate(Prior = Prior^0.1 / sum(Prior^0.1)) %>% 
  mutate(Likelihood = bernoulli_likelihood(theta = theta, data = small_data)) %>% 
  mutate(marginal_likelihood = sum(Prior * Likelihood)) %>% 
  mutate(Posterior = (Prior * Likelihood) / marginal_likelihood) %>% 
  select(theta, Prior, Likelihood, Posterior) %>% 
  gather(key, value, -theta) %>% 
  mutate(key = factor(key, levels = c("Prior", "Likelihood", "Posterior"))) %>% 
  
  ggplot(aes(x = theta, y = value)) +
  geom_area(fill = "grey67") +
  scale_x_continuous(expression(theta), breaks = seq(from = 0, to = 1, by = .2)) +
  ylab("probability density") +
  theme(panel.grid = element_blank()) +
  facet_wrap(~ key, scales = "free_y", ncol = 1)


large_data <- rep(0:1, times = c(30, 10))

tibble(theta = seq(from = 0, to = 1, by = .001),
       Prior = c(seq(from = 0, to = 1, length.out = 501),
                 seq(from = 0.998, to = 0, length.out = 500))) %>% 
  mutate(Prior      = Prior / sum(Prior),
         Likelihood = bernoulli_likelihood(theta = theta, data = large_data)) %>% 
  # here's the important line of code
  mutate(Prior = Prior^10) %>% 
  mutate(marginal_likelihood = sum(Prior * Likelihood)) %>% 
  mutate(Posterior = (Prior * Likelihood) / marginal_likelihood) %>% 
  select(theta, Prior, Likelihood, Posterior) %>% 
  gather(key, value, -theta) %>% 
  mutate(key = factor(key, levels = c("Prior", "Likelihood", "Posterior"))) %>% 
  
  ggplot(aes(x = theta, y = value)) +
  geom_area(fill = "grey67") +
  scale_x_continuous(expression(theta), breaks = seq(from = 0, to = 1, by = .2)) +
  ylab("probability density") +
  theme(panel.grid = element_blank()) +
  facet_wrap(~ key, scales = "free_y", ncol = 1)