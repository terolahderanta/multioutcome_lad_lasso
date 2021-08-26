
library(tidyverse)
library(patchwork)
library(LaplacesDemon)
library(targets)
source("functions.R")

#### RUN TARGETS ####
tar_make()

#### ANALYSIS FUNCTIONS ####

plot_hist_x <- function(dat, title = ""){
  qplot(x = dat$X, 
        geom = "histogram",
        xlab = "Value of the explanatory variable",
        ylab = "Count",
        bins = 20) +
    coord_cartesian(xlim = c(0, 5),
                    ylim = c(0, 420)) +
  ggtitle(title)
   
}

plot_hist_y <- function(dat, title = ""){
  qplot(x = dat$Y[,1], 
        geom = "histogram",
        xlab = "Value of the explanatory variable Y1",
        ylab = "Count",
        bins = 30) +
    #coord_cartesian(xlim = c(-5.5, 10),
    #                ylim = c(0, 22.5)) +
    ggtitle(title)
  
}

table_correct <- function(beta_tibble, prop_zero, error_type, total_correct = 2000){
  beta_tibble |> 
    filter(method != "original") |>
    group_by(method) |> 
    summarise("Total number of correct" = sum(is_correct),
              "Percentage of correct" = sum(is_correct) / total_correct) |>
    rename(Method = method) |> 
    ungroup() |> 
    add_column("Proportion of zeroes" = prop_zero, .after = 1) |> 
    add_column("Type of error" = error_type, .before = 1)
}

#### DATA PLOTS ####

# Explaining variables (does not depend on the error type)
plot_hist_x(tar_read(dat_uniform1), "Proportion of zeros = 0.1") +
plot_hist_x(tar_read(dat_uniform2), "Proportion of zeros = 0.2") +
plot_hist_x(tar_read(dat_uniform3), "Proportion of zeros = 0.3") +
plot_hist_x(tar_read(dat_uniform4), "Proportion of zeros = 0.4")
#ggsave(filename = "images/figdatexm.eps",width = 8, height = 5)

# Response variable
plot_hist_y(tar_read(dat_uniform1), "Proportion of zeros 0.1") +
plot_hist_y(tar_read(dat_uniform3), "Proportion of zeros 0.3") +
plot_hist_y(tar_read(dat_uniform2), "Proportion of zeros 0.2") +
plot_hist_y(tar_read(dat_uniform4), "Proportion of zeros 0.4") +
  plot_annotation(title = "Simulation with uniform error")

plot_hist_y(tar_read(dat_laplace1), "Proportion of zeros 0.1") +
plot_hist_y(tar_read(dat_laplace2), "Proportion of zeros 0.2") +
plot_hist_y(tar_read(dat_laplace3), "Proportion of zeros 0.3") +
plot_hist_y(tar_read(dat_laplace4), "Proportion of zeros 0.4") +
  plot_annotation(title = "Simulation with asymmetric laplace error")


#### RESULT TABLES ####

# Table of how accurate the method was
correct_tibble1 <-
  rbind(
    table_correct(tar_read(beta_uniform1), 0.1, "Uniform"),
    table_correct(tar_read(beta_uniform2), 0.2, "Uniform"),
    table_correct(tar_read(beta_uniform3), 0.3, "Uniform"),
    table_correct(tar_read(beta_uniform4), 0.4, "Uniform"),
    table_correct(tar_read(beta_laplace1), 0.1, "Laplace"),
    table_correct(tar_read(beta_laplace2), 0.2, "Laplace"),
    table_correct(tar_read(beta_laplace3), 0.3, "Laplace"),
    table_correct(tar_read(beta_laplace4), 0.4, "Laplace")
  ) |> 
  add_column(Data = "n > q") 

correct_tibble2 <-
  rbind(
    table_correct(tar_read(beta_uniform1_2), 0.1, "Uniform"),
    table_correct(tar_read(beta_uniform2_2), 0.2, "Uniform"),
    table_correct(tar_read(beta_uniform3_2), 0.3, "Uniform"),
    table_correct(tar_read(beta_uniform4_2), 0.4, "Uniform"),
    table_correct(tar_read(beta_laplace1_2), 0.1, "Laplace"),
    table_correct(tar_read(beta_laplace2_2), 0.2, "Laplace"),
    table_correct(tar_read(beta_laplace3_2), 0.3, "Laplace"),
    table_correct(tar_read(beta_laplace4_2), 0.4, "Laplace")
  ) |> 
  add_column(Data = "q > n")


rbind(correct_tibble1, correct_tibble2) |> 
ggplot(aes(x = `Proportion of zeroes`, 
           y = `Percentage of correct`,
           fill = `Method`)) +
  geom_bar( stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0,1)) +
  scale_y_continuous(position = "right")+
  geom_hline(aes(yintercept = 1)) +
  scale_fill_discrete(labels = c("Adap. LAD-lasso", "LAD-lasso")) +
  facet_grid(cols = vars(`Type of error`), rows = vars(Data), switch = "y") +
  theme(legend.position = "bottom",
        strip.text.y.left = element_text(angle = 0))+
  ylab("Percentage of correct\n")
ggsave(filename = "images/figsim.eps",width = 8, height = 5)

correct_tibble1 |> 
  #filter(`Type of error` == "Uniform") |> 
  ggplot(aes(x = `Proportion of zeroes`, 
             y = `Percentage of correct`,
             fill = `Method`)) +
  geom_bar( stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0,1)) +
  geom_hline(aes(yintercept = 1)) +
  scale_fill_discrete(labels = c("Adap. LAD-lasso", "LAD-lasso")) +
  facet_grid(cols = vars(`Type of error`))

correct_tibble2 |> 
  #filter(`Type of error` == "Uniform") |> 
  ggplot(aes(x = `Proportion of zeroes`, 
             y = `Percentage of correct`,
             fill = `Method`)) +
  geom_bar( stat = "identity", position = position_dodge()) +
  coord_cartesian(ylim = c(0,1))+
  geom_hline(aes(yintercept = 1)) +
  scale_fill_discrete(labels = c("Adap. LAD-lasso", "LAD-lasso"))+
  facet_grid(cols = vars(`Type of error`))




# Plot the estimates on the first two data set
plot_beta_bar(beta_tibble_uniform |> filter(lap == 1))
plot_beta_bar(beta_tibble_uniform |> filter(lap == 2))
plot_beta_bar(beta_tibble_uniform |> filter(lap == 3))
plot_beta_bar(beta_tibble_uniform |> filter(lap == 4))
plot_beta_bar(beta_tibble_uniform |> filter(lap == 5))

 
temp <- fit_n_models(N = 10, 
             error_type = "uniform", 
             n = 100, 
             n_beta = 10, 
             n_effect = 3, 
             n_response = 2, 
             p_zeroes = 0.1)
