######################################## L45

# open dataset

#saveRDS(df_coords, file = "dataset/df_coords.rds")
df_coords <- readRDS("dataset/df_coords.rds")

#saveRDS(db_long, file = "dataset/db_long.rds")
db_long <- readRDS("dataset/db_long.rds")


db_long$sr_name <- substr(db_long$sr_name, 1, 13)

db_long <- db_long %>%
  tidyr::separate(sr_name, into = c("survey", "idnova"), sep = "_")


db_long1 <- db_long %>%
  filter(str_starts(survey, "L45"))


db_long3 <- merge(db_long1,df_coords,by = "ag_name",all.x = T)



db_long3$z <- db_long3$titer
db_long3$variant <- db_long3$ag_name


data <- db_long3
data$z <- as.numeric(data$z)

unique_ids <- unique(data$idnova)
list_fits <- list()

for (id in unique_ids) {
  subdata <- data[data$idnova == id, ]
  fit <- gam(z ~ s(x, y, k = 12), data = subdata)
  list_fits[[id]] <- fit
}



library(mgcv)
library(dplyr)
library(purrr)
library(plotly)

# Datos: data con columnas idnova, x, y, z

# Crear grilla común para predicción
x_seq <- seq(min(data$x), max(data$x), length.out = 50)
y_seq <- seq(min(data$y), max(data$y), length.out = 50)
grid <- expand.grid(x = x_seq, y = y_seq)

# Ajustar GAM para cada idnova
fits <- data %>%
  group_by(idnova) %>%
  group_split() %>%
  set_names(unique(data$idnova)) %>%
  map(~ gam(z ~ s(x, y, k = 12), data = .x))

# Predecir en la grilla para cada modelo
predictions <- map(fits, function(mod) {
  grid$z_hat <- predict(mod, newdata = grid)
  grid
})

# Función para convertir predicciones a matriz para add_surface
make_z_matrix <- function(df) {
  matrix(df$z_hat, nrow = length(x_seq), ncol = length(y_seq))
}

# Graficar cada superficie en plotly, guardar en lista
plots <- map2(predictions, names(predictions), function(df, id) {
  z_mat <- make_z_matrix(df)
  plot_ly() %>%
    add_surface(x = x_seq, y = y_seq, z = z_mat,
                colorscale = "Viridis", showscale = FALSE) %>%
    layout(title = paste("Surface for idnova", id),
           scene = list(
             xaxis = list(title = "x"),
             yaxis = list(title = "y"),
             zaxis = list(title = "z", range = c(0, max(df$z_hat, na.rm = TRUE)))
           ))
})

# Mostrar las gráficas una a una (ejemplo con el primero)
plots[[1]]




# Calcular volumen para cada idnova
volumen_por_idnova <- data %>%
  group_by(idnova) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$idnova))) %>%
  map_dbl(calc_gam_volume_percent)

# Resultado en data frame
volumen_df <- tibble(
  idnova = names(volumen_por_idnova),
  volumen_percent = volumen_por_idnova
)

# Ver
print(volumen_df)


df45 <- readRDS("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/1. Landscape clean/MSD/dataset/volumenes/L45_volumenx2.rds")

df45$volumen100 <- df45$volumen/270*100
df45$id <- substr(df45$id, 5, nchar(df45$id))

df45.2 <- merge(df45,volumen_df,by.x = "id",by.y = "idnova",all = T)


library(ggplot2)
library(ggpubr)
ggplot(df45.2, aes(x = volumen100, y = volumen_percent)) +
  geom_point(color = "darkblue", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  stat_cor(method = "pearson", label.x = min(df45.2$volumen100, na.rm = TRUE),
           label.y = max(df45.2$volumen_percent, na.rm = TRUE), 
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_minimal() +
  labs(
    title = "Correlación entre volumen100 y volumen_percent",
    x = "Volumen 100 (Valor original)",
    y = "Volumen % bajo la superficie GAM"
  )



# Crear figura base
fig <- plot_ly()

# Usar color gris translúcido
color_plomo_trans <- 'rgba(80,80,80,0.3)'

# Añadir todas las superficies (192)
for (i in seq_along(predictions)) {
  z_mat <- make_z_matrix(predictions[[i]])
  
  fig <- fig %>%
    add_surface(
      x = x_seq,
      y = y_seq,
      z = z_mat,
      surfacecolor = matrix(1, nrow = length(x_seq), ncol = length(y_seq)),
      colorscale = list(c(0, 1), c(color_plomo_trans, color_plomo_trans)),
      showscale = FALSE,
      opacity = 0.3,
      name = paste("idnova", names(predictions)[i]),
      showlegend = FALSE
    )
}

# Layout final
fig <- fig %>%
  layout(
    title = "Superficies GAM combinadas (n = 192)",
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z", range = c(0, 1))
    )
  )

fig


######################################## L46


# open dataset

#saveRDS(df_coords, file = "dataset/df_coords.rds")
df_coords <- readRDS("dataset/df_coords.rds")

#saveRDS(db_long, file = "dataset/db_long.rds")
db_long <- readRDS("dataset/db_long.rds")


db_long$sr_name <- substr(db_long$sr_name, 1, 13)

db_long <- db_long %>%
  tidyr::separate(sr_name, into = c("survey", "idnova"), sep = "_")


db_long1 <- db_long %>%
  filter(str_starts(survey, "L46"))


db_long3 <- merge(db_long1,df_coords,by = "ag_name",all.x = T)



db_long3$z <- db_long3$titer
db_long3$variant <- db_long3$ag_name


data <- db_long3
data$z <- as.numeric(data$z)

unique_ids <- unique(data$idnova)
list_fits <- list()

for (id in unique_ids) {
  subdata <- data[data$idnova == id, ]
  fit <- gam(z ~ s(x, y, k = 12), data = subdata)
  list_fits[[id]] <- fit
}



library(mgcv)
library(dplyr)
library(purrr)
library(plotly)

# Datos: data con columnas idnova, x, y, z

# Crear grilla común para predicción
x_seq <- seq(min(data$x), max(data$x), length.out = 50)
y_seq <- seq(min(data$y), max(data$y), length.out = 50)
grid <- expand.grid(x = x_seq, y = y_seq)

# Ajustar GAM para cada idnova
fits <- data %>%
  group_by(idnova) %>%
  group_split() %>%
  set_names(unique(data$idnova)) %>%
  map(~ gam(z ~ s(x, y, k = 12), data = .x))

# Predecir en la grilla para cada modelo
predictions <- map(fits, function(mod) {
  grid$z_hat <- predict(mod, newdata = grid)
  grid
})

# Función para convertir predicciones a matriz para add_surface
make_z_matrix <- function(df) {
  matrix(df$z_hat, nrow = length(x_seq), ncol = length(y_seq))
}

# Graficar cada superficie en plotly, guardar en lista
plots <- map2(predictions, names(predictions), function(df, id) {
  z_mat <- make_z_matrix(df)
  plot_ly() %>%
    add_surface(x = x_seq, y = y_seq, z = z_mat,
                colorscale = "Viridis", showscale = FALSE) %>%
    layout(title = paste("Surface for idnova", id),
           scene = list(
             xaxis = list(title = "x"),
             yaxis = list(title = "y"),
             zaxis = list(title = "z", range = c(0, max(df$z_hat, na.rm = TRUE)))
           ))
})

# Mostrar las gráficas una a una (ejemplo con el primero)
plots[[1]]




# Calcular volumen para cada idnova
volumen_por_idnova <- data %>%
  group_by(idnova) %>%
  group_split() %>%
  set_names(map_chr(., ~ unique(.x$idnova))) %>%
  map_dbl(calc_gam_volume_percent)

# Resultado en data frame
volumen_df <- tibble(
  idnova = names(volumen_por_idnova),
  volumen_percent = volumen_por_idnova
)

# Ver
print(volumen_df)
df46$slope

df46 <- readRDS("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/1. Landscape clean/MSD/dataset/volumenes/L46_volumenx2.rds")

df46$volumen100 <- df46$volumen/270*100
df46$id <- substr(df46$id, 5, nchar(df46$id))

df46.2 <- merge(df46,volumen_df,by.x = "id",by.y = "idnova",all = T)


library(ggplot2)
library(ggpubr)
ggplot(df46.2, aes(x = volumen100, y = volumen_percent)) +
  geom_point(color = "darkblue", size = 2, alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  stat_cor(method = "pearson", label.x = min(df46.2$volumen100, na.rm = TRUE),
           label.y = max(df46.2$volumen_percent, na.rm = TRUE), 
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  theme_minimal() +
  labs(
    title = "Correlación entre volumen100 y volumen_percent",
    x = "Volumen 100 (Valor original)",
    y = "Volumen % bajo la superficie GAM"
  )



# Crear figura base
fig <- plot_ly()

# Usar color gris translúcido
#color_plomo_trans <- 'rgba(80,80,80,0.3)'

# Añadir todas las superficies (192)
for (i in seq_along(predictions)) {
  z_mat <- make_z_matrix(predictions[[i]])
  
  fig <- fig %>%
    add_surface(
      x = x_seq,
      y = y_seq,
      z = z_mat,
      surfacecolor = matrix(1, nrow = length(x_seq), ncol = length(y_seq)),
      colorscale = list(c(0, 1), c(color_plomo_trans, color_plomo_trans)),
      showscale = FALSE,
      opacity = 0.3,
      name = paste("idnova", names(predictions)[i]),
      showlegend = FALSE
    )
}

# Layout final
fig <- fig %>%
  layout(
    title = "Superficies GAM combinadas (n = 192)",
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z", range = c(0, 1))
    )
  )

fig
