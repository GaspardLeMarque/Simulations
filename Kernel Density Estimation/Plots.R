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

IMSE <- rbind(
  data.frame(variable = "standard_normal", value = IMSE_normal, method = "Plug-in bandwidth"),
  data.frame(variable = "standard_normal", value = IMSE_normal_silver, method = "Silvermans rule"),
  data.frame(variable = "standard_normal", value = IMSE_normal_AKDE, method ="Adaptive KDE"),  
  data.frame(variable = "standard_normal", value = IMSE_standard_normal_KNN, method = "kNN"),
  data.frame(variable = "t_distribution", value = IMSE_t_dist, method = "Plug-in bandwidth"),
  data.frame(variable = "t_distribution", value = IMSE_silver_t_dist, method = "Silvermans rule"),
  data.frame(variable = "t_distribution", value = IMSE_t_AKDE, method = "Adaptive KDE"),
  data.frame(variable = "t_distribution", value = IMSE_t_KNN, method = "kNN"),
  data.frame(variable = "bimodal_normal", value = IMSE_bimodal_normal, method = "Plug-in bandwidth"),
  data.frame(variable = "bimodal_normal", value = IMSE_bimodal_normal_silver, method = "Silvermans rule"),
  data.frame(variable = "bimodal_normal", value = IMSE_bimodal_AKDE, method = "Adaptive KDE"),
  data.frame(variable = "bimodal_normal", value = IMSE_bimodal_KNN, method = "kNN"),
  data.frame(variable = "beta_gamma", value = IMSE_gamma_beta, method = "Plug-in bandwidth"),
  data.frame(variable = "beta_gamma", value = IMSE_gamma_beta_silver, method = "Silvermans rule"),
  data.frame(variable = "beta_gamma", value = IMSE_gamma_beta_AKDE, method = "Adaptive KDE"),
  data.frame(variable = "beta_gamma", value = IMSE_gamma_beta_KNN, method = "kNN")
)
IMSE$variable <- factor(x = IMSE$variable, levels = c("standard_normal", "t_distribution", "bimodal_normal", "beta_gamma"))
IMSE$method <- factor(x = IMSE$method, levels = c("Plug-in bandwidth", "Silvermans rule", "Adaptive KDE", "kNN"))

#Boxplot
plot_ly(data = IMSE, x = ~variable, y = ~value, color = ~method, type = "box") %>%
  layout(boxmode = "group", 
         xaxis = list(title="", 
                      showline = FALSE,
                      showticklabels = FALSE,
                      showgrid = FALSE), 
         yaxis = list(title="MISE", titlefont = list(size = 20), tickfont = list(size = 20)),
         annotations = list(x=c(0, 1, 2, 3), 
                            y = c(0.22, 0.22, 0.22, 0.22), 
                            text=sprintf("<b>%s</b>",c("Standard\nNormal", "Bimodal\nNormal", "t-Distribution", "Beta-Gamma", "")),
                            showarrow = FALSE,
                            font = list(size = 20)),
         legend = list(font=list(size=20), x= 0, y=0.85)) #change fontsize (probably 20)

#Benchmark
mean(extract_numeric(log.txt_bimodal_AKDE))
mean(extract_numeric(log.txt_bimodal_KNN))
mean(extract_numeric(log.txt_bimodal))
mean(extract_numeric(log.txt_bimodal_sil))

mean(extract_numeric(log.txt_gamma_beta_AKDE))
mean(extract_numeric(log.txt_gamma_beta_KNN))
mean(extract_numeric(log.txt_gamma_beta))
mean(extract_numeric(log.txt_gamma_beta_sil))

mean(extract_numeric(log.txt_normal_AKDE))
mean(extract_numeric(log.txt_standard_normal_KNN))
mean(extract_numeric(log.txt_normal))
mean(extract_numeric(log.txt_normal_silver))

mean(extract_numeric(log.txt_t_AKDE))
mean(extract_numeric(log.txt_t_KNN))
mean(extract_numeric(log.txt_t))
mean(extract_numeric(log.txt_t_silver))
