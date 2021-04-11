#Plots

#Required packages
library(tidyverse)
library(reshape2)
library(plotly)

#Load in files that were saved previously
load("C:/Users/User 1/Desktop/normal.RDA")
load("C:/Users/User 1/Desktop/normal_silver.RDA")
load("C:/Users/User 1/Desktop/gamma_beta.RDA")
load("C:/Users/User 1/Desktop/gamma_beta_silver.RDA")
load("C:/Users/User 1/Desktop/bimodal.RDA")
load("C:/Users/User 1/Desktop/bimodal_silver.RDA")
load("C:/Users/User 1/Desktop/t_standard.RDA")
load("C:/Users/User 1/Desktop/t_silver.RDA")
load("C:/Users/User 1/Desktop/beta_gamma_akde.RDA")
load("C:/Users/User 1/Desktop/standard_normal_AKDE.RDA")
load("C:/Users/User 1/Desktop/standard_normal_KNN.RDA")
load("C:/Users/User 1/Desktop/gamma_beta_KNN.RDA")
load("C:/Users/User 1/Desktop/t_KNN.RDA")
load("C:/Users/User 1/Desktop/bimodal_KNN.RDA")
load("C:/Users/User 1/Desktop/bimodal_AKDE.RDA")
load("C:/Users/User 1/Desktop/t_AKDE.RDA")