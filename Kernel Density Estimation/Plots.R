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

##### Contour and 3D plots #####

### Bimodal Distribution ###

#Estimation
est_bimodal_normal <- melt(aver_est_bimodal_normal)
names(est_bimodal_normal) <- c("x", "y", "Density")
est_bimodal_AKDE <- melt(aver_est_bimodal_AKDE)
names(est_bimodal_AKDE) <- c("x", "y", "Density")
est_bimodal_KNN <- melt(aver_est_bimodal_KNN)
names(est_bimodal_KNN) <- c("x", "y", "Density")

max(est_bimodal$Density)
max(est_bimodal_AKDE$Density)
max(est_bimodal_KNN$Density)

end = 0.009
size_bimodal <- 0.001
plot_est_bimodal_normal <- plot_ly(est_bimodal_normal,  x = ~x, y = ~y, z= ~Density, type = "contour", showscale = FALSE,
                                   colorscale = 'Jet',
                                   contours = list(
                                     start = 0,
                                     end = 0.01,
                                     size = size_bimodal
                                   )) %>% layout(xaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                                 yaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                                 legend = list(font=list(size = 20)))
plot_est_bimodal_AKDE <- plot_ly(est_bimodal_AKDE,  x = ~x, y = ~y, z= ~Density, type = "contour",showscale = FALSE,
                                 colorscale = 'Jet',
                                 contours = list(
                                   start = 0,
                                   end = 0.01,
                                   size = size_bimodal
                                 )) %>% layout(xaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                               yaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                               legend = list(font=list(size = 20)))

plot_est_bimodal_KNN <- plot_ly(est_bimodal_KNN,  x = ~x, y = ~y, z= ~Density, type = "contour",showscale = FALSE,
                                colorscale = 'Jet',
                                contours = list(
                                  start = 0,
                                  end = 0.01,
                                  size = size_bimodal
                                )) %>% layout(xaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                              yaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                              legend = list(font=list(size = 20)))



p<-subplot(plot_est_bimodal_normal, plot_est_bimodal_AKDE, plot_est_bimodal_KNN)
legend <- plot_ly(est_bimodal_KNN,  x = ~x, y = ~y, z= ~Density, type = "contour",
                  colorscale = 'Jet',
                  contours = list(
                    start = 0,
                    end = 0.01,
                    size = size_bimodal
                  )) %>% layout(xaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                yaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                legend = list(font=list(size = 20)))

#RMSE
rmse_bimodal_normal <- melt(rmse_estimates_bimodal_normal)
names(rmse_bimodal_normal) <- c("x", "y", "Density")
rmse_bimodal_silver <- melt(rmse_estimates_bimodal_normal_silver)
names(rmse_bimodal_silver) <- c("x", "y", "Density")
rmse_bimodal_AKDE <- melt(rmse_estimates_bimodal_AKDE)
names(rmse_bimodal_AKDE) <- c("x", "y", "Density")
rmse_bimodal_KNN <- melt(rmse_estimates_bimodal_KNN)
names(rmse_bimodal_KNN) <- c("x", "y", "Density")

max(rmse_bimodal_normal$Density)
max(rmse_bimodal_silver$Density)
max(rmse_bimodal_AKDE$Density)
max(rmse_bimodal_KNN$Density)
size_bimodal_rmse <- 0.005
plot_rmse_bimodal_normal <- plot_ly(rmse_bimodal_normal,  x = ~x, y = ~y, z= ~Density, type = "contour",showscale = FALSE,
                                    colorscale = 'Jet',
                                    contours = list(
                                      start = 0,
                                      end = 0.05,
                                      size = size_bimodal_rmse
                                    )) %>% layout(xaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                                  yaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                                  legend = list(font=list(size = 20)))



plot_rmse_bimodal_AKDE <- plot_ly(rmse_bimodal_AKDE,  x = ~x, y = ~y, z= ~Density, type = "contour",showscale = FALSE,
                                  colorscale = 'Jet',
                                  contours = list(
                                    start = 0,
                                    end = 0.05,
                                    size = size_bimodal_rmse
                                  )) %>% layout(xaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                                yaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                                legend = list(font=list(size = 20)))


plot_rmse_bimodal_KNN <- plot_ly(rmse_bimodal_KNN,  x = ~x, y = ~y, z= ~Density, type = "contour",showscale = FALSE,
                                 colorscale = 'Jet',
                                 contours = list(
                                   start = 0,
                                   end = 0.05,
                                   size = size_bimodal_rmse
                                 )) %>% layout(xaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                               yaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                               legend = list(font=list(size = 20)))

p<-subplot(plot_rmse_bimodal_normal, plot_rmse_bimodal_AKDE, plot_rmse_bimodal_KNN)
legend <- plot_ly(rmse_bimodal_KNN,  x = ~x, y = ~y, z= ~Density, type = "contour",
                  colorscale = 'Jet',
                  contours = list(
                    start = 0,
                    end = 0.05,
                    size = size_bimodal_rmse
                  )) %>% layout(xaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                yaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                legend = list(font=list(size = 20)))

#Bias
bias_bimodal_normal <- melt(bias_estimates_bimodal_normal)
names(bias_bimodal_normal) <- c("x", "y", "Density")
bias_bimodal_silver <- melt(bias_estimates_bimodal_normal_silver)
names(bias_bimodal_silver) <- c("x", "y", "Density")
bias_bimodal_AKDE <- melt(bias_estimates_bimodal_AKDE)
names(bias_bimodal_AKDE) <- c("x", "y", "Density")
bias_bimodal_KNN <- melt(bias_estimates_bimodal_KNN)
names(bias_bimodal_KNN) <- c("x", "y", "Density")

max(bias_bimodal_normal$Density)
max(bias_bimodal_silver$Density)
max(bias_bimodal_AKDE$Density)
max(bias_bimodal_KNN$Density)

min(bias_bimodal_normal$Density)
min(bias_bimodal_silver$Density)
min(bias_bimodal_AKDE$Density)
min(bias_bimodal_KNN$Density)

start_bias <- -0.046
end_bias <- 0.008
size_bimodal_bias <- 0.01

#Grid resolution for densities
nx = 128
ny = nx
xmin = -10
ymin = -10
xmax = 10
ymax = 10
xcoordinates = seq(xmin, xmax, length.out = nx)
ycoordinates = seq(ymin, ymax, length.out = ny)

z_axis_bimodal <- list(range = c(-0.05,0.008))

plot_bias_bimodal_normal <- plot_ly(
  x = xcoordinates, 
  y = ycoordinates, 
  z = bias_estimates_bimodal_normal, showscale = FALSE, scene = "scene1", colorscale = "Jet",
) %>% 
  add_surface(colorbar=list(title="Density")) %>%
  layout(scene = list(zaxis = z_axis_bimodal))

plot_bias_bimodal_AKDE <- plot_ly(
  x = xcoordinates, 
  y = ycoordinates, 
  z = bias_estimates_bimodal_AKDE, showscale = FALSE,scene = "scene2", colorscale = "Jet",
) %>% 
  add_surface(colorbar=list(title="Density")) %>%
  layout(scene = list(zaxis = z_axis_bimodal))

plot_bias_bimodal_KNN <- plot_ly(
  x = xcoordinates, 
  y = ycoordinates, 
  z = bias_estimates_bimodal_KNN, showscale = FALSE,scene = "scene3", colorscale = "Jet",
) %>% 
  add_surface(colorbar=list(title="Density")) %>%
  layout(scene = list(zaxis = z_axis_bimodal))

p<- subplot(plot_bias_bimodal_normal, plot_bias_bimodal_AKDE, plot_bias_bimodal_KNN)
legend <- plot_ly(
  x = xcoordinates, 
  y = ycoordinates, 
  z = bias_estimates_bimodal_normal,  colorscale = "Jet",
) %>% 
  add_surface(colorbar=list(title="Density")) %>%
  layout(scene = list(zaxis = z_axis_bimodal))

#True density
true_bimodal_normal <- melt(true_est_bimodal_normal)
names(true_bimodal_normal) <- c("x", "y", "Density")
max(true_bimodal_normal$Density)

p <- plot_ly(true_bimodal_normal, x = ~x, y = ~y, z= ~Density, type = "contour", showscale = FALSE,
             colorscale = 'Jet',
             contours = list(
               start = 0,
               end = 0.055,
               size = 0.005
             )) %>% layout(xaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                           yaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                           legend = list(font=list(size = 20)))

legend <- plot_ly(true_bimodal_normal, x = ~x, y = ~y, z= ~Density, type = "contour", 
                  colorscale = 'Jet',
                  contours = list(
                    start = 0,
                    end = 0.055,
                    size = 0.005
                  )) %>% layout(xaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                yaxis = list(tickfont = list(size = 20), titlefont = list(size = 20), showticklabels = FALSE, showgrid = FALSE, showline = FALSE),
                                legend = list(font=list(size = 20)))

p<- plot_ly(
  x = xcoordinates, 
  y = ycoordinates, 
  z = true_est_bimodal_normal, showscale = FALSE, colorscale = "Jet",
) %>% 
  add_surface(colorbar=list(title="Density"))
legend <- plot_ly(
  x = xcoordinates, 
  y = ycoordinates, 
  z = true_est_bimodal_normal, colorscale = "Jet",
) %>% 
  add_surface(colorbar=list(title="Density"))
