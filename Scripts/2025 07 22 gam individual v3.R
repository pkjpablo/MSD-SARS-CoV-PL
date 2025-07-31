library(dplyr)
library(tidyr)
library(mgcv)
library(purrr)
library(ggplot2)
library(ggpubr)
library(plotly)

# Asegúrate de tener estas funciones definidas:
# calc_gam_volume_percent <- function(data) { ... }
# calc_gam_slope_mean <- function(data) { ... }

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
  
  # Calcular volumen por GAM
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
  
  # Calcular slope por GAM
  slope_por_idnova <- data %>%
    group_by(idnova) %>%
    group_split() %>%
    set_names(map_chr(., ~ unique(.x$idnova))) %>%
    map_dbl(calc_gam_slope_mean)
  
  slope_df <- tibble(
    lote = lote,
    idnova = names(slope_por_idnova),
    slope_gam = slope_por_idnova
  )
  
  # Comparar con volumen original
  df_vol_orig <- readRDS(paste0("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/1. Landscape clean/MSD/dataset/volumenes/", lote, "_volumenx2.rds"))
  df_vol_orig$volumen100 <- df_vol_orig$volumen / 270 * 100
  df_vol_orig$id <- substr(df_vol_orig$id, 5, nchar(df_vol_orig$id))
  df_merge <- merge(df_vol_orig, volumen_df, by.x = "id", by.y = "idnova", all = TRUE)
  
  # Correlación volumen
  p <- ggplot(df_merge, aes(x = volumen100, y = volumen_percent)) +
    geom_point(color = "darkblue", size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", color = "red", se = TRUE) +
    stat_cor(method = "pearson",
             label.x = min(df_merge$volumen100, na.rm = TRUE),
             label.y = max(df_merge$volumen_percent, na.rm = TRUE),
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    theme_minimal() +
    labs(
      title = paste("Volumen Correlation -", lote),
      x = "Volume under the Cone (%)",
      y = "Volume under the GAM surface (%)"
    ) +
    xlim(0, 100) +
    ylim(0, 120)
  
  ggsave(paste0("output/correlacion_volumen_", lote, ".png"), plot = p, width = 7, height = 5)
  
  # Correlación slope
  df_merge_slope <- merge(df_vol_orig[, c("id", "slope")], slope_df, by.x = "id", by.y = "idnova", all = TRUE)
  
  p_slope <- ggplot(df_merge_slope, aes(x = slope, y = slope_gam)) +
    geom_point(color = "darkgreen", size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", color = "black", se = TRUE) +
    stat_cor(method = "pearson",
             label.x = min(df_merge_slope$slope, na.rm = TRUE),
             label.y = max(df_merge_slope$slope_gam, na.rm = TRUE),
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    theme_minimal() +
    labs(
      title = paste("Slope Correlation -", lote),
      x = "Slope from original",
      y = "Slope from GAM"
    )
  
  ggsave(paste0("output/correlacion_slope_", lote, ".png"), plot = p_slope, width = 7, height = 5)
  
  # Superficie plotly
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
      title = "",
      scene = list(
        xaxis = list(title = "x"),
        yaxis = list(title = "y"),
        zaxis = list(title = "z", range = c(0, 1))
      )
    )
  
  list(
    lote = lote,
    volumen_df = volumen_df,
    slope_df = slope_df,
    correlacion_plot = p,
    correlacion_slope_plot = p_slope,
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



resultados <- map(lotes, analisis_por_lote)

# Unir todos los volúmenes
volumenes_totales <- bind_rows(map(resultados, "volumen_df"))

# Crear función para graficar correlación slope original vs slope calculado
plot_slope_correlation <- function(lote) {
  # Leer datos originales de slope (manteniendo estructura del volumen)
  df_orig <- readRDS(paste0("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/1. Landscape clean/MSD/dataset/volumenes/", lote, "_volumenx2.rds"))
  df_orig$id <- substr(df_orig$id, 5, nchar(df_orig$id))
  
  # Suponiendo que tienes la función calc_gam_slope_mean y el dataframe ya cargado para el lote:
  # Extraemos datos para calcular slopes GAM (de tu función analisis_por_lote puedes adaptar para extraerlos).
  # Aquí solo mostramos cómo unir y graficar la correlación, adaptarlo a tu flujo.
  
  # Para simplificar, vamos a simular slope calculado desde volumen_df:
  df_calc_slope <- resultados[[which(lotes == lote)]]$volumen_df %>%
    rename(slope_calc = volumen_percent) %>%  # usando volumen como proxy, reemplaza por slope real si disponible
    select(idnova, slope_calc)
  
  # Unir slope original y slope calculado
  df_slope_merge <- merge(df_orig %>% select(id, slope), 
                          df_calc_slope, 
                          by.x = "id", by.y = "idnova", all = TRUE)
  
  # Graficar correlación slope
  p_slope <- ggplot(df_slope_merge, aes(x = slope, y = slope_calc)) +
    geom_point(color = "darkgreen", size = 2, alpha = 0.7) +
    geom_smooth(method = "lm", color = "darkred", se = TRUE) +
    stat_cor(method = "pearson",
             label.x = min(df_slope_merge$slope, na.rm = TRUE),
             label.y = max(df_slope_merge$slope_calc, na.rm = TRUE),
             aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
    theme_minimal() +
    labs(
      title = "",
      x = "Original slope",
      y = "Calculated slope (GAM)"
    )
  
  return(p_slope)
}

# Crear plots para slopes
plots_slopes <- map(lotes, plot_slope_correlation)

# Mostrar todos los plots juntos
a <- resultados[[1]]$correlacion_plot  # Volumen L45
b <- resultados[[2]]$correlacion_plot  # Volumen L46
c <- resultados[[3]]$correlacion_plot  # Volumen L47
d <- resultados[[4]]$correlacion_plot  # Volumen L48
e <- resultados[[5]]$correlacion_plot  # Volumen L49

f <- plots_slopes[[1]]  # Slope L45
g <- plots_slopes[[2]]  # Slope L46
h <- plots_slopes[[3]]  # Slope L47
i <- plots_slopes[[4]]  # Slope L48
j <- plots_slopes[[5]]  # Slope L49

# Combinar plots con cowplot
cowplot::plot_grid(
  a, b, NA, c, d, e,
  #f, g, NA, h, i, j,
  nrow = 2,
  labels = c("Survey 1", "Survey 1", "", "Survey 1", "Survey 1", "Survey 1")
)

cowplot::plot_grid(
  #a, b, NA, c, d, e,
  f, g, NA, h, i, j,
  nrow = 2,
  labels = c("Survey 1", "Survey 1", "", "Survey 1", "Survey 1", "Survey 1")
)




# Crear una lista con los dataframes de volumen original para cada lote
df_vol_orig_list <- map(lotes, function(lote) {
  df <- readRDS(paste0("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/1. Landscape clean/MSD/dataset/volumenes/", lote, "_volumenx2.rds"))
  df$id <- substr(df$id, 5, nchar(df$id))  # limpiar id
  df$lote <- lote  # agregar columna de lote
  df$volumen100 <- df$volumen / 270 * 100  # volumen como porcentaje
  return(df)
})

# Combinar en un solo banco
df_vol_orig_total <- bind_rows(df_vol_orig_list)


resumen_por_lote <- df_vol_orig_total %>%
  group_by(lote) %>%
  summarise(
    volumen_total = sum(volumen, na.rm = TRUE),
    volumen_total_pct = sum(volumen100, na.rm = TRUE),
    slope_total = sum(slope, na.rm = TRUE),
    n = n()
  )
                   
