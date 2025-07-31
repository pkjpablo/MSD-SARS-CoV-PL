#save.acmap(full_map_p1_adj, "maps/PL p23 and 32 L45-L48 all v2.ace")
setwd("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/1. Landscape clean/MSD")
# ------------------------------------------ landscape --------------------------------------

### Use the version 4.2.2 

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

### colores del mapa
mapColors <- read.csv("dataset/Ejemplo 3/map-colors v2.csv", 
                      comment = ";",row.names = 'Antigen', header = TRUE)

# read the full map roessler 

map_orig =  read.acmap("maps/pau da lima p23 and 32 L45-49 all.ace")
plot(map_orig)

### colores de los grupos
sr_colors <- read_delim("dataset/Ejemplo 3/sr_group_colors.csv", 
                        delim = ";", 
                        escape_double = FALSE, 
                        trim_ws = TRUE)
sr_colors <- sr_colors %>%
  tibble::column_to_rownames(var = "SerumGroup")

# read the full map PL to obtain the titers 

map_jp <-  read.acmap("maps/pau da lima p23 and 32 L45-49 all.ace")
plot(map_jp)
map_long <- long_map_info(map_jp)


# Selec a survey
#map_long <-map_long[grep("^L45", map_long$sr_name), ]
#map_long <-map_long[grep("^L46", map_long$sr_name), ]
#map_long <-map_long[grep("^L49", map_long$sr_name), ]
#map_long <-map_long[grep("^L48", map_long$sr_name), ]
map_long <-map_long[grep("^L49", map_long$sr_name), ]

head(map_long)
table(map_long$sr_name)

# Using 0.01 for negative values p23 for ancestral and p32 for BA.1
hist(as.numeric(map_long$titer))

library(dplyr)

## exponencial de base 2
map_long <- map_long %>%
  mutate(titer = 10*2^(as.numeric(titer)*10))
head(map_long)
map_long %>% ggplot(aes(x=titer))+geom_histogram()+scale_x_continuous(trans = "log")

map_long %>%
  select(titer, ag_name, sr_name, sr_group) -> titerdata

titerdata


#titerdata$sr_group <- "L45"
#titerdata$sr_group <- substr(titerdata$sr_name, 1, 3)
titerdata$sr_group <- substr(titerdata$sr_name, 1,13 )
#titerdata$sr_group <- titerdata$sr_name


dbtiter <- titerdata %>%
  pivot_wider(
    names_from = ag_name,  # Las columnas en el formato ancho se crearán a partir de ag_name
    values_from = titer     # Los valores de esas columnas vendrán de la columna titer
  )


grupos <- split(dbtiter[, c("sr_name", setdiff(names(dbtiter), c("sr_group", "sr_name")))], dbtiter$sr_group)

matrices_list <- lapply(grupos, function(df) {
  as.matrix(df)
})

titertables <- lapply(matrices_list, function(mat) {
  rownames(mat) <- mat[, 1]  # Asignar sr_name como nombres de fila
  mat <- mat[, -1, drop = FALSE]  # Eliminar la columna de sr_name
  return(mat)
})


lndscp_fits <- lapply(
  titertables,
  function(titertable) {
    
    ablandscape.fit(
      titers = titertable, #[,ags_to_fit_lndscp],
      bandwidth = 1,
      degree = 1,
      method = "cone",
      error.sd = 1,
      acmap = full_map_p1_adj,
      control = list(
        optimise.cone.slope = TRUE
      )
    )
  }
)


#titertables_groups <- group_data(titerdata)




#titertables_groups <- group_data(titerdata[,c("sr_group","ag_name")])
titertables_groups <- group_by(titerdata,sr_group)

grouped_tibble <- titerdata[,c("sr_group","ag_name")] %>%
  group_by(sr_group) %>%
  dplyr::summarise(.rows = list(row_number()), .groups = "drop")
grouped_tibble

grouped_positions <- titerdata %>%
  group_by(sr_group) %>%
  dplyr::summarize(.rows = list(which(sr_group == first(sr_group))), .groups = "drop")

result <- titerdata %>%
  group_by(sr_group) %>%
  dplyr::summarise(.rows = list(which(row_number() %in% row_number())), .groups = "drop")

# Verificar si el archivo existe
if (file.exists("dataset/L49x_gmt.rds")) {
  # Cargar gmt_data del archivo
  gmt_data <- readRDS(file = "dataset/L49x_gmt.rds")
} else {
  # Calcular gmt_data y guardarlo si el archivo no existe
  gmt_data <- titerdata %>%
    group_by(
      sr_group,
      ag_name
    ) %>%
    dplyr::summarize(gmt = titertools::gmt(titer, dilution_stepsize = 0)["mean", "estimate"])
  
  # Guardar gmt_data en el archivo
  saveRDS(gmt_data, file = "dataset/L49x_gmt.rds")
}

# Mostrar gmt_data
print(gmt_data)

# angle for html page
angle <- list(
  rotation = c(-1.4427, 0.0100, -0.0263), #c(-1.3365, 0.0055, -0.0576),# c(-1.4592, 0.0045, -0.0144)
  translation = c(0, 0,0), #translation = c(0.0344, 0.0459, 0.1175),
  zoom = 1.5
  # zoom = 1.1646 # higher is more zoomed out
)


#Modified function for base plot
base_plot_data3js_grs <- function(map, lndscp_fits, highlighted_ags, lims, ag_plot_names, alternative_ba5 = FALSE, opti_nr = 1,
                                  add_border = TRUE, add_axis = TRUE){
  
  x_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 1])
  y_coords <- c(agCoords(map)[agNames(map) %in% highlighted_ags, 2])
  z_coords <- rep(0.02, length(highlighted_ags))
  ag_point_size <- c(rep(14, length(highlighted_ags))) / 5
  ag_col <- c(agOutline(map)[agNames(map) %in% highlighted_ags])
  ag_fill <- c(agFill(map)[agNames(map) %in% highlighted_ags])
  labels <- c(ag_plot_names[agNames(map) %in% highlighted_ags])
  border_col <- "grey50"
  
  z_lims <- c(0,10)
  axis_at <- seq(z_lims[1], z_lims[2],2)
  # Setup plot
  data3js <- ablandscapes:::lndscp3d_setup(
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = z_lims,
    aspect.z = 0.5,
    options = list(
      lwd.grid =  0.05,
      sidegrid.lwd = 1,
      sidegrid.col = border_col,
      sidegrid.at = list("z" = axis_at),
      zaxt = "lin"
    ),
    show.axis = FALSE
  )
  
  if(add_axis){
    
    axis_labels <- seq(0,1,by=0.2)
    
    data3js <- r3js::axis3js(
      data3js,
      side = "z",
      at = axis_at,
      labels = axis_labels,
      # labeloffset = 0.11,
      cornerside = "f",
      size = 20,
      alignment = "right"
    )
  }
  
  # Add basemap
  data3js <- lndscp3d_map(
    data3js = data3js,
    fit = lndscp_fits[[1]],
    xlim = lims$xlim,
    ylim = lims$ylim,
    zlim = c(0, 10),
    show.map.sera = FALSE,
    options = list(
      opacity.basemap = 0.3
    )
  )
  
  data3js <- r3js::points3js(
    data3js,
    x          = x_coords,
    y          = y_coords,
    z          = z_coords,
    size       = ag_point_size,
    col        = ag_col,
    fill       = ag_fill,
    lwd        = 0.5,
    opacity    = 1,
    highlight  = list(col = "red"),
    label      = labels,
    toggle     = "Basepoints",
    depthWrite = FALSE,
    shape      = "circle filled"
  )
  
  if(add_border){
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[1]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[2],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    # y border
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[1], lims$ylim[1]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    data3js <- lines3js(data3js, x = c(lims$xlim[1],lims$xlim[2]), y = c(lims$ylim[2], lims$ylim[2]), z = c(0, 0),
                        lwd = 1.2, col = border_col)
    
    data3js <- r3js::box3js(
      data3js,
      col   = border_col
    )
    
  }
  
  return(data3js)
}

gam_fit$residuals
gam_fit$aic
summary(gam_fit)

summary(lndscp_fits)

#Base plot
lndscp_list <- list()
data3js <- base_plot_data3js_grs(map, lndscp_fits, agNames(map), lims = list(xlim=c(-2,6),ylim=c(-3,2)), agNames(map), add_axis = T)
data3js



# do plot 
PlotL49 <- plot_landscapes_from_list(data3js,
                                     titertables_groups,
                                     lndscp_fits,
                                     map,
                                     gmt_data,
                                     agNames(map),
                                     agNames(map),
                                     lndscp_colors = sr_colors,
                                     show_gmts = F)



# Verificar si el archivo existe
if (file.exists("dataset/PlotL49x.rds")) {
  # Cargar gmt_data del archivo
  gmt_data <- readRDS(file = "dataset/PlotL49x.rds")
} else {
  # Crear gmt_data utilizando la función si el archivo no existe
  PlotL49 <- plot_landscapes_from_list(data3js, titertables_groups, lndscp_fits, map, gmt_data, 
                                       agNames(map), agNames(map), lndscp_colors = sr_colors,
                                       show_gmts = FALSE)
  # Guardar gmt_data en un archivo .rds para futuros usos
  saveRDS(PlotL49, file = "dataset/PlotL49x.rds")
  message("El archivo 'PlotL49x.rds' no se encontró. Se ha creado y guardado un nuevo archivo.")
}


PlotL49


######## Calcualr volumenes e slope

# Inicializar una lista para almacenar los resultados
rm(df)
df <- data.frame(id = character(), slope = numeric(), stringsAsFactors = FALSE)

# Iterar sobre los elementos de lndscp_fits
for (i in seq_along(lndscp_fits)) {
  # Obtener el ID del elemento actual
  id <- names(lndscp_fits)[i]
  
  # Obtener el slope del elemento actual
  slope <- lndscp_fits[[i]][["cone"]][["cone_slope"]]
  
  # Agregar una nueva fila al data frame
  df <- rbind(df, data.frame(id = id, slope = slope))
}

print(df)




# Extraer índices únicos (tomando uno de cada par) para obetner z

result <- sapply(PlotL49[["plot"]], function(elem) {
  !is.null(elem$x) && !is.null(elem$y) && !is.null(elem$z)
})
unique_indices <- which(result)[seq(1, length(which(result)), by = 2)]


result <- sapply(PlotL49[["plot"]], function(elem) {
  !is.null(elem$x) && !is.null(elem$y) && !is.null(elem$z)
})
unique_indices <- which(result)[seq(1, length(which(result)), by = 2)]



# Definir la función para calcular el volumen, considerando z negativo como 0
calcular_volumen <- function(index) {
  x <- PlotL49[["plot"]][[unique_indices[index]]][["x"]]
  y <- PlotL49[["plot"]][[unique_indices[index]]][["y"]]
  z <- PlotL49[["plot"]][[unique_indices[index]]][["z"]]
  
  # Establecer a 0 todos los valores negativos de z
  z[z < 0] <- 0
  
  # Convertir x, y, z a matrices de 13x10
  x_matrix <- matrix(x, nrow = 10, ncol = 13, byrow = TRUE)
  y_matrix <- matrix(y, nrow = 10, ncol = 13, byrow = TRUE)
  z_matrix <- matrix(z, nrow = 10, ncol = 13, byrow = TRUE)
  
  # Calcular el área de cada celda
  dx <- abs(x_matrix[1, 1] - x_matrix[1, 2])  # Distancia en x
  dy <- abs(y_matrix[1, 1] - y_matrix[2, 1])  # Distancia en y
  
  # Calcular el volumen usando la regla del trapecio
  volumen <- 0
  
  for (i in 1:(nrow(z_matrix) - 1)) {
    for (j in 1:(ncol(z_matrix) - 1)) {
      # Calcular el volumen de cada celda
      volume_cell <- (z_matrix[i, j] + z_matrix[i + 1, j] + z_matrix[i, j + 1] + z_matrix[i + 1, j + 1]) / 4 * dx * dy
      volumen <- volumen + volume_cell
    }
  }
  
  return(volumen)
}




# Aplicar la función a todos los índices en unique_indices, con valor absoluto si se desea
volumenes <- sapply(1:length(unique_indices), function(i) calcular_volumen(i))

# Imprimir todos los volúmenes
print(volumenes)

# Calcular el porcentaje
porcentaje_volumenes <- volumenes / 270 * 100

# Calcular la media y el error estándar
media <- mean(porcentaje_volumenes)
error_estandar <- sd(porcentaje_volumenes) / sqrt(length(porcentaje_volumenes))

# Calcular el intervalo de confianza del 95%
IC_inf <- media - 1.96 * error_estandar
IC_sup <- media + 1.96 * error_estandar

# Imprimir resultados
cat("Media:", media, "\n")
cat("IC 95%:", IC_inf, "-", IC_sup, "\n")




# Agregar los volúmenes al dataframe
df$volumen <- volumenes
dfL49 <- df

# 

# Inicializar una lista vacía para almacenar los valores extraídos
z_values_list <- list()

# Recorrer cada índice en unique_indices
for (i in unique_indices) {
  # Extraer los cuatro valores de z en las posiciones especificadas
  z_1_1 <- PlotL49[["plot"]][[i]][["z"]][1, 1]
  z_1_10 <- PlotL49[["plot"]][[i]][["z"]][1, 10]
  z_13_1 <- PlotL49[["plot"]][[i]][["z"]][13, 1]
  z_13_10 <- PlotL49[["plot"]][[i]][["z"]][13, 10]
  
  # Almacenar los valores en la lista
  z_values_list[[i]] <- c(z_1_1, z_1_10, z_13_1, z_13_10)
}

# Convertir la lista a un data.frame
z_values_df <- as.data.frame(do.call(rbind, z_values_list))

# Opcional: Asignar nombres a las columnas
colnames(z_values_df) <- c("z_1_1", "z_1_10", "z_13_1", "z_13_10")

# Mostrar el data.frame
print(z_values_df)


if (nrow(dfL49) == nrow(z_values_df)) {
  # Unir ambos data.frames por las columnas
  dfL49 <- cbind(dfL49, z_values_df)
  
  # Mostrar el data.frame combinado
  print(dfL49)
} else {
  print("Los data.frames no tienen el mismo número de filas.")
}


# Definir la ruta del archivo
file_path <- "dataset/volumenes/L49_volumenx2.rds"

# Comprobar si el archivo ya existe en la ubicación indicada
if (!file.exists(file_path)) {
  # Si el archivo no existe, guardar el archivo
  saveRDS(dfL49, file = file_path)
  cat("El archivo se ha guardado correctamente.\n")
} else {
  # Si el archivo existe, no hacer nada
  cat("El archivo ya existe y no se ha guardado.\n")
}





# Gráfico para z_1_1, escalado y con rango ajustado (iz)
ggplot(dfL49, aes(x = factor(1), y = z_1_1/10)) + 
  geom_violin(fill = "skyblue", color = "darkblue") +
  labs(title = "", x = "", y = "") +
  coord_cartesian(ylim = c(-0.5, 1)) +  # Limitar el rango del eje Y
  theme_minimal()

# Gráfico para z_1_10, escalado y con rango ajustado (atras)
ggplot(dfL49, aes(x = factor(1), y = z_1_10/10)) + 
  geom_violin(fill = "skyblue", color = "darkblue") +
  labs(title = "", x = "", y = "") +
  coord_cartesian(ylim = c(-0.5, 1)) +  # Limitar el rango del eje Y
  theme_minimal()+ theme_classic()

# Gráfico para z_13_1, escalado y con rango ajustado (frente)
ggplot(dfL49, aes(x = factor(1), y = z_13_1/10)) + 
  geom_violin(fill = "skyblue", color = "darkblue") +
  labs(title = "", x = "", y = "") +
  coord_cartesian(ylim = c(-0.5, 1)) +  # Limitar el rango del eje Y
  theme_minimal()+ theme_classic()

# Gráfico para z_13_10, escalado y con rango ajustado (derecha)
ggplot(dfL49, aes(x = factor(1), y = z_13_10/10)) + 
  geom_violin(fill = "skyblue", color = "darkblue") +
  labs(title = "", x = "", y = "") +
  coord_cartesian(ylim = c(-0.5, 1)) +  # Limitar el rango del eje Y
  theme_minimal()+ theme_classic()

261/270


# Cargar la librería plotly para gráficos interactivos 3D
library(plotly)

# Crear un gráfico tridimensional utilizando plotly
fig <- plot_ly(dfL49, 
               x = ~z_1_1/10, y = ~1, z = ~1,
               color = ~1, colors = c('skyblue', 'darkblue'), 
               type = "violin", mode = "markers") %>%
  layout(title = "Gráfico Tridimensional de las Variables Z",
         scene = list(
           xaxis = list(title = "z_1_1/10"),
           yaxis = list(title = "z_1_10/10"),
           zaxis = list(title = "z_13_1/10"),
           coloraxis = list(colorbar = list(title = "z_13_10/10"))
         ))

# Mostrar el gráfico
fig





# Lista para almacenar los índices de objetos donde z es negativo
negative_z_indices <- list()

# Iterar sobre cada índice en unique_indices para detectar valores negativos en z
for (index in unique_indices) {
  z <- PlotL49[["plot"]][[index]][["z"]]
  
  if (any(z < 0)) {
    negative_z_indices <- c(negative_z_indices, index)
  }
}

# Imprimir los índices de los objetos donde z es negativo
print(negative_z_indices)



PlotL49


# Volumen total


PlotL49[["plot"]][[238]][["x"]]
PlotL49[["plot"]][[238]][["y"]]
PlotL49[["plot"]][[11460]][["x"]]
PlotL49[["plot"]][[11460]][["y"]]

#valor maximo 
(1.09473995 +4.90526005)* (2.913616+1.586384)*10

