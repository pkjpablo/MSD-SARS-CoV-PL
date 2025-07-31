library(dplyr)
library(tidyr)
library(mgcv)
library(purrr)
library(ggplot2)
library(ggpubr)
library(plotly)

# Asegúrate de tener esta función definida
# calc_gam_volume_percent <- function(data) { ... }

# Función general para procesar un lote
analisis_por_lote <- function(lote) {
  message("Procesando ", lote, "...")
  
  # Leer datos
  df_coords <- readRDS("dataset/df_coords.rds")
  db_long <- readRDS("dataset/db_long.rds")
  
  # Preprocesamiento
  db_long$sr_name <- substr(db_long$sr_name, 1, 13)
  db_long <- db_long %>% separate(sr_name, into = c("survey", "idnova"), sep = "_")
  db_long1 <- db_long %>% filter(str_starts(survey, lote))
  db_long3 <- merge(db_long1, df_coords, by = "ag_name", all.x = TRUE)
  
  db_long3$z <- as.numeric(db_long3$titer)
  db_long3$variant <- db_long3$ag_name
  data <- db_long3
  
  # Ajustar GAM
  fits <- data %>%
    group_by(idnova) %>%
    group_split() %>%
    set_names(unique(data$idnova)) %>%
    map(~ gam(z ~ s(x, y, k = 12), data = .x))
  
  # Grilla común
  x_seq <- seq(min(data$x), max(data$x), length.out = 50)
  y_seq <- seq(min(data$y), max(data$y), length.out = 50)
  grid <- expand.grid(x = x_seq, y = y_seq)
  
  # Predicciones
  predictions <- map(fits, function(mod) {
    grid$z_hat <- predict(mod, newdata = grid)
    grid
  })
  
  # Función auxiliar
  make_z_matrix <- function(df) {
    matrix(df$z_hat, nrow = length(x_seq), ncol = length(y_seq))
  }
  
  # Volumen por GAM
  volumen_por_idnova <- data %>%
    group_by(idnova) %>%
    group_split() %>%
    set_names(map_chr(., ~ unique(.x$idnova))) %>%
    map_dbl(calc_gam_volume_percent)
  
  volumen_df <- tibble(
    lote = lote,
    idnova = names(volumen_por_idnova),
    volumen_percent = volumen_por_idnova
  )
  
  # Comparar con volumen original
  df_vol_orig <- readRDS(paste0("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/1. Landscape clean/MSD/dataset/volumenes/", lote, "_volumenx2.rds"))
  df_vol_orig$volumen100 <- df_vol_orig$volumen / 270 * 100
  df_vol_orig$id <- substr(df_vol_orig$id, 5, nchar(df_vol_orig$id))
  df_merge <- merge(df_vol_orig, volumen_df, by.x = "id", by.y = "idnova", all = TRUE)
  
  # Gráfico de correlación
  p <- ggplot(df_merge, aes(x = volumen100, y = volumen_percent)) +
    geom_point(color = "darkblue", size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    stat_cor(method = "pearson",
             label.x = min(df_merge$volumen100, na.rm = TRUE),
             label.y = max(df_merge$volumen_percent, na.rm = TRUE),
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    theme_minimal() +
    labs(
      title = paste("", lote),
      x = "Volume under the Cone (%)",
      y = "Volume under the GAM surface (%)"
    ) +
    xlim(0, 100) +   # Escala fija para eje X
    ylim(0, 120)     # Escala fija para eje Y
  
  # Guardar PNG (opcional)
  ggsave(paste0("output/correlacion_", lote, ".png"), plot = p, width = 7, height = 5)
  
  # Crear figura plotly con todas las superficies
  fig <- plot_ly()
  for (i in seq_along(predictions)) {
    z_mat <- make_z_matrix(predictions[[i]])
    fig <- fig %>%
      add_surface(
        x = x_seq,
        y = y_seq,
        z = z_mat,
        surfacecolor = matrix(1, nrow = length(x_seq), ncol = length(y_seq)),
        colorscale = list(c(0, 1), c("rgba(80,80,80,0.3)", "rgba(80,80,80,0.3)")),
        showscale = FALSE,
        opacity = 0.3,
        name = paste("idnova", names(predictions)[i]),
        showlegend = FALSE
      )
  }
  
  fig <- fig %>%
    layout(
      title = paste("Superficies GAM combinadas -", lote),
      scene = list(
        xaxis = list(title = "x"),
        yaxis = list(title = "y"),
        zaxis = list(title = "z", range = c(0, 1))
      )
    )
  
  list(
    lote = lote,
    volumen_df = volumen_df,
    correlacion_plot = p,
    superficie_plotly = fig
  )
}

# Lotes a procesar
lotes <- c("L45", "L46", "L47", "L48", "L49")

# Ejecutar análisis para todos
resultados <- map(lotes, analisis_por_lote)

# Unir todos los volúmenes
volumenes_totales <- bind_rows(map(resultados, "volumen_df"))


a <- print(resultados[[1]]$correlacion_plot)  # L45
b <- print(resultados[[2]]$correlacion_plot)  # L46
c <- print(resultados[[3]]$correlacion_plot)  # L47
d <- print(resultados[[4]]$correlacion_plot)  # L48
e <- print(resultados[[5]]$correlacion_plot)  # L49

cowplot::plot_grid(a,b,NA,c,d,e, nrow = 2)


walk(resultados, ~ print(.x$correlacion_plot))



resultados[[1]]$superficie_plotly  # L45
resultados[[2]]$superficie_plotly  # L46
resultados[[3]]$superficie_plotly  # L47
resultados[[4]]$superficie_plotly  # L48
resultados[[5]]$superficie_plotly  # L49
