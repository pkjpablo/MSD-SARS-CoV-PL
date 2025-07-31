# open dataset

#saveRDS(df_coords, file = "dataset/df_coords.rds")
df_coords <- readRDS("dataset/df_coords.rds")

#saveRDS(db_long, file = "dataset/db_long.rds")
db_long <- readRDS("dataset/db_long.rds")



# L45

db_long1 <- db_long %>%
  filter(str_starts(sr_name, "L45"))

db_long3 <- merge(db_long1,df_coords,by = "ag_name",all.x = T)



x <- db_long3$x
y <- db_long3$y
z <- db_long3$titer
variant <- db_long3$ag_name


## Dataset simple

data <- data.frame(x = x, y = y, z = z, variant=variant)

str(data$z)
str(data$y)
str(data$x)


data$z <- as.numeric(data$z)



## Base map 

data2 <- data %>%
  group_by(x, y,variant) %>%
  summarise(z = mean(z, na.rm = TRUE), .groups = "drop")


# Varaibles
data2$x
data2$y
data2$variant

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

library(plotly)


data2$variant <- as.character(data2$variant)
data2$color <- mapColors[data2$variant]
data2$y_label <- data2$y + 0.05

fig <- plot_ly() %>%
  # Puntos grandes con color por variante
  add_markers(
    data = data2,
    x = ~x,
    y = ~y,
    color = ~variant,
    colors = mapColors,
    marker = list(size = 12),
    showlegend = FALSE,
    hoverinfo = "text",
    text = ~variant
  ) %>%
  # Etiquetas con el nombre de la variante, justo encima de cada punto
  add_text(
    data = data2,
    x = ~x,
    y = ~y_label,
    text = ~variant,
    textfont = list(size = 12, color = "black"),
    showlegend = FALSE
  ) %>%
  layout(
    title = "Puntos por Variante con Nombres",
    xaxis = list(title = "x"),
    yaxis = list(title = "y")
  )

fig



fig <- plot_ly(
  data = data2,
  x = ~x,
  y = ~y,
  z = ~z,
  type = 'scatter3d',
  mode = 'markers+text',
  marker = list(
    size = 6,             # tamaño de la "esfera"
    #symbol = 'circle',    # símbolo de esfera
    #color = ~variant,     # mapea color a variante
    #colors = "black",   # usa tu paleta personalizada
    opacity = 0.8
  ),
  text = ~variant,         # muestra el nombre de la variante
  textposition = 'top center',
  hoverinfo = 'text'
) %>%
  layout(
    title = "",
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    )
  )

fig



# Crear figura
fig <- plot_ly()

# Agregar puntos planos (proyección en z = 0)
fig <- fig %>%
  add_markers(
    data = data2,
    x = ~x,
    y = ~y,
    z = 0,  # Fijamos z = 0
    color = ~variant,
    colors = mapColors,
    marker = list(size = 8, opacity = 0.3, symbol = "circle"),
    hoverinfo = "text",
    text = ~paste("Base -", variant),
    name = "Base (z = 0)"
  )

# Agregar puntos reales en 3D
fig <- fig %>%
  add_markers(
    data = data2,
    x = ~x,
    y = ~y,
    z = ~z,
    color = ~variant,
    colors = mapColors,
    marker = list(size = 6, opacity = 0.9),
    hoverinfo = "text",
    text = ~variant,
    name = "Puntos 3D"
  )

# Agregar etiquetas sobre los puntos 3D
fig <- fig %>%
  add_text(
    data = data2,
    x = ~x,
    y = ~y,
    z = ~z,
    text = ~variant,
    textposition = 'top center',
    textfont = list(size = 10, color = "black"),
    showlegend = FALSE
  )

# Agregar líneas negras entre base y punto 3D
for (i in 1:nrow(data2)) {
  fig <- fig %>%
    add_trace(
      x = c(data2$x[i], data2$x[i]),
      y = c(data2$y[i], data2$y[i]),
      z = c(0, data2$z[i]),
      type = "scatter3d",
      mode = "lines",
      line = list(color = "black", width = 2),
      showlegend = FALSE
    )
}

# Layout final
fig <- fig %>%
  layout(
    title = "Puntos en 3D con proyección y líneas de conexión",
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    )
  )

fig


## gam 


# --- PREPARA DATOS ---

data3 <- data

# Ajuste GAM para superficie
gam_fit <- gam(z ~ s(x, y, k = 12), data = data3)

x_seq <- seq(min(data3$x), max(data3$x), length.out = 100)
y_seq <- seq(min(data3$y), max(data3$y), length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)
grid$z_hat <- predict(gam_fit, newdata = grid)

z_matrix <- matrix(grid$z_hat, nrow = 100, ncol = 100)
z_matrix <- t(z_matrix)[ , ncol(z_matrix):1]
z_matrix <- z_matrix[nrow(z_matrix):1, ncol(z_matrix):1]

# --- INICIA FIGURA COMBINADA ---
fig <- plot_ly()

# Superficie GAM
fig <- fig %>%
  add_surface(
    x = x_seq,
    y = y_seq,
    z = z_matrix,
    opacity = 0.6,
    showscale = FALSE,
    colorscale = "Viridis",
    name = "Superficie GAM"
  )

# Puntos planos base (z = 0)
fig <- fig %>%
  add_markers(
    data = data2,
    x = ~x,
    y = ~y,
    z = 0,
    color = ~variant,
    colors = mapColors,
    marker = list(size = 8, opacity = 0.3, symbol = "circle"),
    hoverinfo = "text",
    text = ~paste("Base -", variant),
    name = "Base (z = 0)"
  )

# Puntos reales 3D
fig <- fig %>%
  add_markers(
    data = data2,
    x = ~x,
    y = ~y,
    z = ~z,
    color = ~variant,
    colors = mapColors,
    marker = list(size = 6, opacity = 0.9),
    hoverinfo = "text",
    text = ~variant,
    name = "Puntos 3D"
  )

# Etiquetas de los puntos
fig <- fig %>%
  add_text(
    data = data2,
    x = ~x,
    y = ~y,
    z = ~z,
    text = ~variant,
    textposition = 'top center',
    textfont = list(size = 10, color = "black"),
    showlegend = FALSE
  )

# Líneas negras entre puntos base y puntos 3D
for (i in 1:nrow(data2)) {
  fig <- fig %>%
    add_trace(
      x = c(data2$x[i], data2$x[i]),
      y = c(data2$y[i], data2$y[i]),
      z = c(0, data2$z[i]),
      type = "scatter3d",
      mode = "lines",
      line = list(color = "black", width = 2),
      showlegend = FALSE
    )
}

# Layout final
fig <- fig %>%
  layout(
    title = "Survey 1",
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z")
    )
  )

fig

fig <- fig %>%
  layout(
    title = "Survey 1",
    scene = list(
      xaxis = list(title = "x"),
      yaxis = list(title = "y"),
      zaxis = list(title = "z", range = c(0, 1))  # <- aquí fijamos el rango
    )
  )


fig













data <- data2

# Fit 2D thin plate spline
gam_fit <- gam(z ~ s(x, y, k = 12), data = data)

# Create grid for predictions
x_seq <- seq(min(data$x), max(data$x), length.out = 100)
y_seq <- seq(min(data$y), max(data$y), length.out = 100)
grid <- expand.grid(x = x_seq, y = y_seq)

grid$z_hat <- predict(gam_fit, newdata = grid)




# Reshape prediction surface into matrix form for 3D plotting
z_matrix <- matrix(grid$z_hat, nrow = 100, ncol = 100)

z_matrix <- t(z_matrix)[, nrow(z_matrix):1]
#z_matrix <- z_matrix[nrow(z_matrix):1, ncol(z_matrix):1]
#z_matrix <- t(z_matrix[nrow(z_matrix):1, ])



# Ajuste en las coordenadas z para colocar el texto encima de los puntos
data$z_text <- data$z + 0.1  # Ajusta 0.1 a tu preferencia para mover el texto hacia arriba

# Gráfico 3D con plotly
L45p <- plot_ly(x = x_seq, y = y_seq, z = ~z_matrix) %>%
  add_surface(opacity = 0.8) %>%
  add_markers(data = data, 
              x = ~x, y = ~y, z = ~z,
              marker = list(color = 'red', size = 5),
              name = "Observed Points",
              text = ~variant,  # Muestra el nombre de la variante
              textposition = 'middle center',  # Asegura que el texto esté centrado encima del marcador
              showlegend = FALSE) %>%
  add_text(data = data, 
           x = ~x, y = ~y, z = ~z_text,  # Usamos z_text para posicionar el texto sobre los puntos
           text = ~variant,  # Agregar el texto de la variante
           showlegend = FALSE, 
           textfont = list(size = 10, color = "blue")) %>%  # Personaliza el texto
  layout(title = "3D Surface Plot of GAM Fit with Observed Points",
         scene = list(xaxis = list(title = "x"),
                      yaxis = list(title = "y"),
                      zaxis = list(title = "Observed & Fitted z")))

L45p
