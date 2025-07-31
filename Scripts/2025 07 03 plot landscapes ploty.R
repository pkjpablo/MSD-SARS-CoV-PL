#setup page and load metadata
#rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
library(titertools)
library(ablandscapes) 
library(r3js)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)

set.seed(100)

source("./functions/remove_reactivity_bias.R")
source("./functions/map_longinfo.R")
source("./functions/sams_landscape_functions.R")


agCoords(full_map_p1_adj)

df_coords <- as.data.frame(agCoords(full_map_p1_adj))

df_coords$ag_name <- rownames(df_coords)     # Convertir los nombres de fila en una columna
rownames(df_coords) <- NULL               # Remover los nombres de fila
colnames(df_coords)[1:2] <- c("x", "y")   # Renombrar las columnas



agCoords(full_map_p1_adj)

plot(full_map_p1_adj)



map_jp <- read.acmap("maps/pau da lima p23 and 32 L45-49 all.ace")
plot(map_jp)
map_long <- long_map_info(map_jp)


db_long <- data.frame(
  titer = map_long$titer,
  sr_name = map_long$sr_name,
  ag_name = map_long$ag_name
)

db_long1 <- db_long %>%
  filter(str_starts(sr_name, "L45"))



db_long3 <- merge(db_long1,df_coords,by = "ag_name",all.x = T)


x <- db_long3$x
y <- db_long3$y
z <- db_long3$titer
data <- data.frame(x = x, y = y, z = z)



library(mgcv)
library(plotly)
library(dplyr)



data$z <- as.numeric(data$z)

# Fit 2D thin plate spline
gam_fit <- gam(z ~ s(x, y, k = 12), data = data)

# Create grid for predictions
x_seq <- seq(min(data$x), max(data$x), length.out = 100)
y_seq <- seq(min(data$y), max(data$y), length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)
grid$z_hat <- predict(gam_fit, newdata = grid)

# Reshape prediction surface into matrix form for 3D plotting
z_matrix <- matrix(grid$z_hat, nrow = 100, ncol = 100)

# Plot with plotly
L45p <-  plot_ly(x = x_seq, y = y_seq, z = ~z_matrix) %>%
  add_surface() %>%
  layout(title = "3D Surface Plot of GAM Fit",
         scene = list(xaxis = list(title = "x"),
                      yaxis = list(title = "y"),
                      zaxis = list(title = "Predicted z")))

L45p

nome <- db_long3$ag_name    # nombre de cada punto

set.seed(123)  # para reproducibilidad
df_plot <- data.frame(
  x = jitter(x, amount = 0.3),
  y = jitter(y, amount = 0.3),
  z = z,
  nome = nome
)


df_labels <- df_plot %>%
  group_by(nome) %>%
  filter(z == max(z, na.rm = TRUE)) %>%
  slice(1) %>%  # en caso de empate
  ungroup()


fig <- plot_ly(
  data = df_plot,
  x = ~x, y = ~y, z = ~z,
  type = "scatter3d",
  mode = "markers",
  marker = list(
    size = 4,
    color = ~z,
    colorscale = "Viridis",
    showscale = TRUE
  )
)

fig <- fig %>%
  add_trace(
    data = df_labels,
    x = ~x, y = ~y, z = ~z,
    type = "scatter3d",
    mode = "text",
    text = ~nome,
    textposition = "top center",
    textfont = list(color = "black", size = 12),
    showlegend = FALSE
  ) %>%
  layout(
    scene = list(
      xaxis = list(title = "y (db_long3$y)"),
      yaxis = list(title = "x (db_long3$x)"),
      zaxis = list(title = "titer")
    )
  )

fig



db_long1 <- db_long %>%
  filter(str_starts(sr_name, "L46"))



db_long3 <- merge(db_long1,df_coords,by = "ag_name",all.x = T)


x <- db_long3$x
y <- db_long3$y
z <- db_long3$titer
data <- data.frame(x = x, y = y, z = z)



library(mgcv)
library(plotly)
library(dplyr)


str(data$x)
str(data$y)
str(data$z)

data$z <- as.numeric(data$z)

# Fit 2D thin plate spline
gam_fit <- gam(z ~ s(x, y, k = 12), data = data)

# Create grid for predictions
x_seq <- seq(min(data$x), max(data$x), length.out = 100)
y_seq <- seq(min(data$y), max(data$y), length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)
grid$z_hat <- predict(gam_fit, newdata = grid)

# Reshape prediction surface into matrix form for 3D plotting
z_matrix <- matrix(grid$z_hat, nrow = 100, ncol = 100)

# Plot with plotly
L46p <-  plot_ly(x = x_seq, y = y_seq, z = ~z_matrix) %>%
  add_surface() %>%
  layout(title = "3D Surface Plot of GAM Fit",
         scene = list(xaxis = list(title = "x"),
                      yaxis = list(title = "y"),
                      zaxis = list(title = "Predicted z")))

L46p




db_long1 <- db_long %>%
  filter(str_starts(sr_name, "L47"))



db_long3 <- merge(db_long1,df_coords,by = "ag_name",all.x = T)


x <- db_long3$x
y <- db_long3$y
z <- db_long3$titer
data <- data.frame(x = x, y = y, z = z)



library(mgcv)
library(plotly)
library(dplyr)


str(data$x)
str(data$y)
str(data$z)

data$z <- as.numeric(data$z)

# Fit 2D thin plate spline
gam_fit <- gam(z ~ s(x, y, k = 12), data = data)

# Create grid for predictions
x_seq <- seq(min(data$x), max(data$x), length.out = 100)
y_seq <- seq(min(data$y), max(data$y), length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)
grid$z_hat <- predict(gam_fit, newdata = grid)

# Reshape prediction surface into matrix form for 3D plotting
z_matrix <- matrix(grid$z_hat, nrow = 100, ncol = 100)

# Plot with plotly
L47p <-  plot_ly(x = x_seq, y = y_seq, z = ~z_matrix) %>%
  add_surface() %>%
  layout(title = "3D Surface Plot of GAM Fit",
         scene = list(xaxis = list(title = "x"),
                      yaxis = list(title = "y"),
                      zaxis = list(title = "Predicted z")))


L47p



db_long1 <- db_long %>%
  filter(str_starts(sr_name, "L48"))



db_long3 <- merge(db_long1,df_coords,by = "ag_name",all.x = T)


x <- db_long3$x
y <- db_long3$y
z <- db_long3$titer
data <- data.frame(x = x, y = y, z = z)



library(mgcv)
library(plotly)
library(dplyr)


str(data$x)
str(data$y)
str(data$z)

data$z <- as.numeric(data$z)

# Fit 2D thin plate spline
gam_fit <- gam(z ~ s(x, y, k = 12), data = data)

# Create grid for predictions
x_seq <- seq(min(data$x), max(data$x), length.out = 100)
y_seq <- seq(min(data$y), max(data$y), length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)
grid$z_hat <- predict(gam_fit, newdata = grid)

# Reshape prediction surface into matrix form for 3D plotting
z_matrix <- matrix(grid$z_hat, nrow = 100, ncol = 100)

# Plot with plotly
L48p <-  plot_ly(x = x_seq, y = y_seq, z = ~z_matrix) %>%
  add_surface() %>%
  layout(title = "3D Surface Plot of GAM Fit",
         scene = list(xaxis = list(title = "x"),
                      yaxis = list(title = "y"),
                      zaxis = list(title = "Predicted z")))

L48p


db_long1 <- db_long %>%
  filter(str_starts(sr_name, "L49"))


db_long3 <- merge(db_long1,df_coords,by = "ag_name",all.x = T)


x <- db_long3$x
y <- db_long3$y
z <- db_long3$titer
data <- data.frame(x = x, y = y, z = z)



library(mgcv)
library(plotly)
library(dplyr)


str(data$x)
str(data$y)
str(data$z)

data$z <- as.numeric(data$z)

# Fit 2D thin plate spline
gam_fit <- gam(z ~ s(x, y, k = 12), data = data)

# Create grid for predictions
x_seq <- seq(min(data$x), max(data$x), length.out = 100)
y_seq <- seq(min(data$y), max(data$y), length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)
grid$z_hat <- predict(gam_fit, newdata = grid)

# Reshape prediction surface into matrix form for 3D plotting
z_matrix <- matrix(grid$z_hat, nrow = 100, ncol = 100)

# Plot with plotly
L49p <-  plot_ly(x = x_seq, y = y_seq, z = ~z_matrix) %>%
  add_surface() %>%
  layout(title = "3D Surface Plot of GAM Fit",
         scene = list(xaxis = list(title = "x"),
                      yaxis = list(title = "y"),
                      zaxis = list(title = "Predicted z")))

L45p
L46p
L47p
L48p
L49p

AIC(gam_fit)


summary(gam_fit)



