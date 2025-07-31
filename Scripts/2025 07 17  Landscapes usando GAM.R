# Pablo 
# SARS-CoV-2 Landscape GAM
# surveys 


## Libraries ----

library(dplyr)
library(mgcv)
library(plotly)
library(stringr)



## Dataset----
#Roessler map
df_coords <- readRDS("dataset/df_coords.rds")

# pl dataset
db_long <- readRDS("dataset/db_long.rds")


## Colores de las variantes ----

mapColors <- c(
  "D614G" = "#393b79",
  "BF.7" = "#647A39",
  "B.1.617.2" = "#d18652",
  "XBB.1" = "#5B004C",
  "BA.2" = "#5B004C",
  "B.1.1.7" = "#637939",
  "BQ.1.3" = "#CD9B1D",
  "P.1.1" = "#b471ab",
  "BA.1" = "#EF3737",
  "BA.5.2.1" = "#bfdaa0",
  "B.1.351" = "#107f97",
  "BQ.1.1" = "#107f97"
)




df_coords <- df_coords %>%
  mutate(ag_name = recode(ag_name,
                          "D614G"      = "Ancestral",
                          "B.1.351"    = "Beta",
                          "P.1.1"      = "Gamma",
                          "B.1.1.7"    = "Alpha",
                          "B.1.617.2"  = "Delta",
                          "BA.5.2.1"   = "BA.5",
                          "BQ.1.3"     = "BQ.1",
                          "BQ.1.1"     = "BQ.1.1",
                          "BF.7"       = "BF.7",
                          "XBB.1"      = "XBB.1",
                          "BA.1"       = "BA.1",
                          "BA.2"       = "BA.2"))




db_long <- db_long %>%
  mutate(ag_name = recode(ag_name,
                          "D614G"      = "Ancestral",
                          "B.1.351"    = "Beta",
                          "P.1.1"      = "Gamma",
                          "B.1.1.7"    = "Alpha",
                          "B.1.617.2"  = "Delta",
                          "BA.5.2.1"   = "BA.5",
                          "BQ.1.3"     = "BQ.1"))

mapColors <- c(
  "Ancestral" = "#393b79",
  "Alpha"     = "#637939",  # B.1.1.7
  "Beta"      = "#CD9B1D",  # B.1.351
  "Gamma"     = "#b471ab",  # P.1.1
  "Delta"     = "#d18652",  # B.1.617.2
  "BA.1"      = "#EF3737",
  "BA.2"      = "#5B004C",
  "BA.5"      = "#bfdaa0",  # BA.5.2.1
  "BF.7"      = "#647A39",
  "BQ.1"      = "#107f97",  # BQ.1.3
  "BQ.1.1"    = "#107f97",
  "XBB.1"     = "#5B004C"
)





## Funcion para crear el landscape con GAM ----

make_survey_fig <- function(
    survey_name, db_long, df_coords,
    z_range = c(0, 100), text_size = 14,
    expand_xy = 0.5, show_labels = TRUE,
    show_ci = TRUE
) {
  db_sub <- db_long %>% filter(str_starts(sr_name, survey_name))
  merged <- merge(db_sub, df_coords, by = "ag_name", all.x = TRUE)
  
  data <- data.frame(
    x = merged$x,
    y = merged$y,
    z = as.numeric(merged$titer),
    variant = as.character(merged$ag_name)
  )
  
  data2 <- data %>%
    group_by(x, y, variant) %>%
    summarise(
      z_mean = mean(z, na.rm = TRUE),
      z_sd = sd(z, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(
      se = z_sd / sqrt(n),
      z_low = z_mean - 1.96 * se,
      z_up = z_mean + 1.96 * se,
      hover_label = sprintf(
        "%s: %.1f (%.1f - %.1f)",
        variant, z_mean, z_low, z_up
      )
    )
  
  gam_fit <- gam(z ~ s(x, y, k = 12), data = data)
  
  x_seq <- seq(min(data$x), max(data$x), length.out = 100)
  y_seq <- seq(min(data$y), max(data$y), length.out = 100)
  grid <- expand.grid(x = x_seq, y = y_seq)
  preds <- predict(gam_fit, newdata = grid, se.fit = TRUE)
  
  grid$z_hat <- preds$fit
  grid$z_low <- preds$fit - 1.96 * preds$se.fit
  grid$z_up <- preds$fit + 1.96 * preds$se.fit
  
  z_matrix <- matrix(grid$z_hat, nrow = length(x_seq), ncol = length(y_seq))
  z_matrix <- t(z_matrix)
  
  z_low_matrix <- matrix(grid$z_low, nrow = length(x_seq), ncol = length(y_seq))
  z_low_matrix <- t(z_low_matrix)
  
  z_up_matrix <- matrix(grid$z_up, nrow = length(x_seq), ncol = length(y_seq))
  z_up_matrix <- t(z_up_matrix)
  
  x_min_exp <- min(data$x) - expand_xy
  x_max_exp <- max(data$x) + expand_xy
  y_min_exp <- min(data$y) - expand_xy
  y_max_exp <- max(data$y) + expand_xy
  
  fig <- plot_ly()
  
  # Superficie principal
  fig <- fig %>%
    add_surface(
      x = x_seq, y = y_seq, z = z_matrix,
      opacity = 0.6, showscale = FALSE, colorscale = "Viridis"
    )
  
  # Superficies de IC (opcional)
  if (show_ci) {
    fig <- fig %>%
      add_surface(
        x = x_seq, y = y_seq, z = z_low_matrix,
        opacity = 0.25, showscale = FALSE, colors = "#999999", name = "IC inferior"
      ) %>%
      add_surface(
        x = x_seq, y = y_seq, z = z_up_matrix,
        opacity = 0.25, showscale = FALSE, colors = "#999999", name = "IC superior"
      )
  }
  
  # Puntos base
  fig <- fig %>%
    add_markers(
      data = data2,
      x = ~x, y = ~y, z = 0,
      color = ~variant,
      colors = mapColors,
      marker = list(size = 8, opacity = 0.3, symbol = "circle"),
      hoverinfo = "text",
      text = ~paste("Base -", variant),
      name = "Base (z = 0)"
    )
  
  # Puntos 3D con media e IC en tooltip
  fig <- fig %>%
    add_markers(
      data = data2,
      x = ~x, y = ~y, z = ~z_mean,
      color = ~variant,
      colors = mapColors,
      marker = list(size = 6, opacity = 0.9),
      hoverinfo = "text",
      text = ~hover_label,
      name = "Pontos 3D"
    )
  
  # Etiquetas (opcional)
  if (show_labels) {
    fig <- fig %>%
      add_text(
        data = data2,
        x = ~x, y = ~y, z = ~z_mean,
        text = ~variant,
        textposition = "top center",
        textfont = list(size = text_size, color = "black"),
        showlegend = FALSE
      )
  }
  
  # Líneas desde el plano base
  for (i in seq_len(nrow(data2))) {
    fig <- fig %>%
      add_trace(
        x = c(data2$x[i], data2$x[i]),
        y = c(data2$y[i], data2$y[i]),
        z = c(0, data2$z_mean[i]),
        type = "scatter3d",
        mode = "lines",
        line = list(color = "black", width = 2),
        showlegend = FALSE
      )
  }
  
  fig <- fig %>%
    layout(
      title = "",
      scene = list(
        xaxis = list(title = "", tickfont = list(size = text_size), range = c(x_min_exp, x_max_exp)),
        yaxis = list(title = "", tickfont = list(size = text_size), range = c(y_min_exp, y_max_exp)),
        zaxis = list(title = "", range = c(0, 1), tickformat = ".0%", tickfont = list(size = text_size)),
        camera = list(eye = list(x = 1.5, y = 1.5, z = 0.7))  # rotación original
      ),
      showlegend = FALSE
    )
  
  return(fig)
}


fig_L45 <- make_survey_fig("L45", db_long, df_coords, show_labels = F, show_ci = TRUE)
fig_L46 <- make_survey_fig("L46", db_long, df_coords, show_labels = F, show_ci = TRUE)
fig_L47 <- make_survey_fig("L47", db_long, df_coords, show_labels = F, show_ci = TRUE)
fig_L48 <- make_survey_fig("L48", db_long, df_coords, show_labels = F, show_ci = TRUE)
fig_L49 <- make_survey_fig("L49", db_long, df_coords, show_labels = F, show_ci = TRUE)


ajustar_figura_3D <- function(figura, zoom = 1.7, x_eye = 0.1, y_eye = -1.25, z_eye = 0.2, text_size = 14) {
  figura %>%
    layout(
      scene = list(
        zaxis = list(
          title = "",
          range = c(0, 1),
          tickmode = "array",
          tickvals = seq(0, 1, 0.2),
          ticktext = paste0(seq(0, 100, 20), "%"),
          tickfont = list(size = text_size)
        ),
        camera = list(
          eye = list(
            x = x_eye * zoom,
            y = y_eye * zoom,
            z = z_eye * zoom
          )
        ),
        xaxis = list(tickfont = list(size = text_size)),
        yaxis = list(tickfont = list(size = text_size))
      )
    )
}



fig_L45v <- ajustar_figura_3D(fig_L45)
fig_L46v <- ajustar_figura_3D(fig_L46)
fig_L47v <- ajustar_figura_3D(fig_L47)
fig_L48v <- ajustar_figura_3D(fig_L48)
fig_L49v <- ajustar_figura_3D(fig_L49)

library(htmltools)

# Crear espacio vacío
empty_div <- tags$div(style = "width:100%; height:400px;")

# Crear panel HTML con estilo CSS para 2x3 layout
panel_layout <- tags$div(
  style = "
    display: grid;
    grid-template-columns: repeat(3, 33%);
    grid-template-rows: auto auto;
    gap: 10px;
  ",
  # Primera fila (2 gráficos + espacio vacío)
  tags$div(fig_L45v), 
  tags$div(fig_L46v),
  empty_div,
  
  # Segunda fila (3 gráficos)
  tags$div(fig_L47v),
  tags$div(fig_L48v),
  tags$div(fig_L49v)
)

# Mostrar como página HTML
browsable(panel_layout)


save_html(browsable(panel_layout), file = "index.html")


## Función para extraer banco de datos por survey
get_survey_data <- function(survey_name, db_long, df_coords) {
  db_sub <- db_long %>% filter(str_starts(sr_name, survey_name))
  merged <- merge(db_sub, df_coords, by = "ag_name", all.x = TRUE)
  
  data <- data.frame(
    variant = merged$ag_name,
    x = merged$x,
    y = merged$y,
    z = as.numeric(merged$titer)
  )
  
  return(data)
}

# Crear bancos de datos para cada survey
df_L45 <- get_survey_data("L45", db_long, df_coords)
df_L46 <- get_survey_data("L46", db_long, df_coords)
df_L47 <- get_survey_data("L47", db_long, df_coords)
df_L48 <- get_survey_data("L48", db_long, df_coords)
df_L49 <- get_survey_data("L49", db_long, df_coords)



#### Volumen 

calc_gam_volume_percent <- function(data) {
  # Ajustar modelo GAM
  gam_fit <- gam(z ~ s(x, y, k = 12), data = data)
  
  # Crear grilla de predicción
  x_seq <- seq(min(data$x), max(data$x), length.out = 100)
  y_seq <- seq(min(data$y), max(data$y), length.out = 100)
  grid <- expand.grid(x = x_seq, y = y_seq)
  grid$z_hat <- predict(gam_fit, newdata = grid)
  
  # Calcular delta de área por celda
  dx <- diff(range(x_seq)) / (length(x_seq) - 1)
  dy <- diff(range(y_seq)) / (length(y_seq) - 1)
  dA <- dx * dy
  
  # Volumen real: suma de z_hat * área de cada celda
  V_real <- sum(grid$z_hat, na.rm = TRUE) * dA
  
  # Volumen máximo: área total * altura máxima (z = 1)
  V_max <- (diff(range(x_seq)) * diff(range(y_seq))) * 1
  
  # Porcentaje
  perc <- (V_real / V_max) * 100
  
  return(perc)
}



vol_L45 <- calc_gam_volume_percent(df_L45)
vol_L46 <- calc_gam_volume_percent(df_L46)
vol_L47 <- calc_gam_volume_percent(df_L47)
vol_L48 <- calc_gam_volume_percent(df_L48)
vol_L49 <- calc_gam_volume_percent(df_L49)


round(c(L45 = vol_L45,
        L46 = vol_L46,
        L47 = vol_L47,
        L48 = vol_L48,
        L49 = vol_L49), 2)


## 95% IC 

calc_gam_volume_percent_ci <- function(data, n_boot = 500) {
  # Almacenar resultados
  boot_volumes <- numeric(n_boot)
  
  for (i in 1:n_boot) {
    # Re-muestreo con reemplazo
    boot_data <- data[sample(1:nrow(data), replace = TRUE), ]
    
    # Ajustar modelo GAM
    gam_fit <- gam(z ~ s(x, y, k = 12), data = boot_data)
    
    # Crear grilla de predicción
    x_seq <- seq(min(data$x), max(data$x), length.out = 100)
    y_seq <- seq(min(data$y), max(data$y), length.out = 100)
    grid <- expand.grid(x = x_seq, y = y_seq)
    grid$z_hat <- predict(gam_fit, newdata = grid)
    
    # Calcular delta de área por celda
    dx <- diff(range(x_seq)) / (length(x_seq) - 1)
    dy <- diff(range(y_seq)) / (length(y_seq) - 1)
    dA <- dx * dy
    
    # Volumen real
    V_real <- sum(grid$z_hat, na.rm = TRUE) * dA
    
    # Volumen máximo
    V_max <- (diff(range(x_seq)) * diff(range(y_seq))) * 1
    
    # Porcentaje
    boot_volumes[i] <- (V_real / V_max) * 100
  }
  
  # Resultado final
  volume_mean <- mean(boot_volumes)
  volume_ci <- quantile(boot_volumes, probs = c(0.025, 0.975))
  
  return(list(
    mean = volume_mean,
    lower = volume_ci[1],
    upper = volume_ci[2]
  ))
}

vol_L45 <- calc_gam_volume_percent_ci(df_L45)
vol_L46 <- calc_gam_volume_percent_ci(df_L46)
vol_L47 <- calc_gam_volume_percent_ci(df_L47)
vol_L48 <- calc_gam_volume_percent_ci(df_L48)
vol_L49 <- calc_gam_volume_percent_ci(df_L49)

vol_L45
vol_L46
vol_L47
vol_L48
vol_L49


library(reactable)

df <- data.frame(
  Lote = c("L45", "L46", "L47", "L48", "L49"),
  Valor = c(
    "8.6  (8.0 – 9.2)",
    "25.9  (24.1 – 27.8)",
    "69.4  (67.5 – 71.2)",
    "76.2  (74.5 – 77.7)",
    "81.7  (80.1 – 83.3)"
  )
)

reactable(df, bordered = TRUE, highlight = TRUE, striped = TRUE)




## pendiente usando gam 

calc_gam_slope_mean <- function(data) {
  # Ajustar modelo GAM
  gam_fit <- gam(z ~ s(x, y, k = 12), data = data)
  
  # Crear grilla
  x_seq <- seq(min(data$x), max(data$x), length.out = 100)
  y_seq <- seq(min(data$y), max(data$y), length.out = 100)
  grid <- expand.grid(x = x_seq, y = y_seq)
  grid$z_hat <- predict(gam_fit, newdata = grid)
  
  # Convertir a matriz
  z_mat <- matrix(grid$z_hat, nrow = length(x_seq), ncol = length(y_seq))
  
  # Calcular diferencias finitas para derivadas en x e y
  dz_dx <- apply(z_mat, 2, diff) / diff(x_seq)[1]
  dz_dy <- t(apply(z_mat, 1, diff)) / diff(y_seq)[1]
  
  # Hacer que tengan misma dimensión (99x100 y 100x99 → 99x99)
  dz_dx <- dz_dx[, -ncol(dz_dx), drop = FALSE]
  dz_dy <- dz_dy[-nrow(dz_dy), , drop = FALSE]
  
  # Calcular magnitud de gradiente (pendiente)
  slope <- sqrt(dz_dx^2 + dz_dy^2)
  
  # Promedio de la pendiente
  mean_slope <- mean(slope, na.rm = TRUE)
  
  return(mean_slope)
}



slope_L45 <- calc_gam_slope_mean(df_L45)
slope_L46 <- calc_gam_slope_mean(df_L46)
slope_L47 <- calc_gam_slope_mean(df_L47)
slope_L48 <- calc_gam_slope_mean(df_L48)
slope_L49 <- calc_gam_slope_mean(df_L49)

round(c(L45 = slope_L45,
        L46 = slope_L46,
        L47 = slope_L47,
        L48 = slope_L48,
        L49 = slope_L49), 4)



# Función para bootstrap del slope promedio
bootstrap_gam_slope <- function(data, n_boot = 500, conf_level = 0.95) {
  slopes <- numeric(n_boot)
  
  for (i in 1:n_boot) {
    sample_data <- data[sample(1:nrow(data), replace = TRUE), ]
    slopes[i] <- tryCatch(calc_gam_slope_mean(sample_data), error = function(e) NA)
  }
  
  # Remover NA por errores de ajuste
  slopes <- slopes[!is.na(slopes)]
  
  # Calcular IC percentil
  alpha <- (1 - conf_level) / 2
  ci_lower <- quantile(slopes, probs = alpha)
  ci_upper <- quantile(slopes, probs = 1 - alpha)
  slope_mean <- mean(slopes)
  
  return(list(
    mean = slope_mean,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    slopes = slopes  # puedes quitar esto si no quieres guardar los valores
  ))
}


boot_L45 <- bootstrap_gam_slope(df_L45)
boot_L46 <- bootstrap_gam_slope(df_L46)
boot_L47 <- bootstrap_gam_slope(df_L47)
boot_L48 <- bootstrap_gam_slope(df_L48)
boot_L49 <- bootstrap_gam_slope(df_L49)

slopes_summary <- data.frame(
  Variante = c("L45", "L46", "L47", "L48", "L49"),
  Slope_Mean = round(c(boot_L45$mean, boot_L46$mean, boot_L47$mean, boot_L48$mean, boot_L49$mean), 4),
  CI_Lower = round(c(boot_L45$ci_lower, boot_L46$ci_lower, boot_L47$ci_lower, boot_L48$ci_lower, boot_L49$ci_lower), 4),
  CI_Upper = round(c(boot_L45$ci_upper, boot_L46$ci_upper, boot_L47$ci_upper, boot_L48$ci_upper, boot_L49$ci_upper), 4)
)

print(slopes_summary)


library(ggplot2)
library(dplyr)

# Datos principales
data <- data.frame(
  date = as.Date(c("2021-01-07", "2021-08-26", "2022-06-12", "2023-02-11", "2024-01-08")),
  volume_mean = c(8.60, 25.86, 69.40, 76.15, 81.7),
  volume_lower = c(8.07, 24.14, 67.55, 74.50, 80.1),
  volume_upper = c(9.20, 27.73, 71.23, 77.82, 83.3),
  Slope_Mean = c(0.0465, 0.0479, 0.0689, 0.0587, 0.0440),
  CI_Lower = c(0.0381, 0.0284, 0.0526, 0.0427, 0.0305),
  CI_Upper = c(0.0546, 0.0725, 0.0884, 0.0757, 0.0577)
)

# Escalar Slope y sus IC para el eje secundario
data <- data %>%
  mutate(
    slope_scaled = Slope_Mean * 500,
    ci_lower_scaled = CI_Lower * 500,
    ci_upper_scaled = CI_Upper * 500
  )

# Fechas de inicio y fin de surveys (pares consecutivos)
vertical_lines <- as.Date(c(
  "2019-09-09", "2019-11-09",
  "2020-11-18", "2021-02-26",
  "2021-07-14", "2021-10-09",
  "2022-03-25", "2022-08-31",
  "2022-11-16", "2023-05-09",
  "2023-10-25", "2024-03-24"
))

# Crear pares de fechas para representar áreas de survey
survey_periods <- data.frame(
  start = vertical_lines[seq(1, length(vertical_lines), by = 2)],
  end = vertical_lines[seq(2, length(vertical_lines), by = 2)]
)

# Gráfico
Fig3 <- ggplot() +
  # Fondo de surveys
  geom_rect(data = survey_periods,
            aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
            fill = "gray90", alpha = 0.4) +
  
  # Área de confianza del Slope
  geom_ribbon(data = data,
              aes(x = date, ymin = ci_lower_scaled, ymax = ci_upper_scaled),
              fill = "#F7B7A3", alpha = 0.3) +
  
  # Área de confianza del Volume
  geom_ribbon(data = data,
              aes(x = date, ymin = volume_lower, ymax = volume_upper),
              fill = "#A3C4BC", alpha = 0.3) +
  
  # Línea y puntos de Volume
  geom_line(data = data,
            aes(x = date, y = volume_mean),
            color = "#A3C4BC", size = 1) +
  geom_point(data = data,
             aes(x = date, y = volume_mean),
             color = "#A3C4BC", size = 3) +
  
  # Línea y puntos de Slope (escalado ×500)
  geom_line(data = data,
            aes(x = date, y = slope_scaled),
            color = "#F7B7A3", size = 1) +
  geom_point(data = data,
             aes(x = date, y = slope_scaled),
             color = "#F7B7A3", size = 3) +
  
  # Ejes
  scale_x_date(date_breaks = "3 months", date_labels = "%b-%Y",
               limits = as.Date(c('2020-03-01', '2024-06-01'))) +
  scale_y_continuous(
    name = "Volume (%)",
    sec.axis = sec_axis(~ . / 500, name = "Slope")
  ) +
  
  # Estilo
  labs(x = "", title = "") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )

# Mostrar gráfico
Fig3






