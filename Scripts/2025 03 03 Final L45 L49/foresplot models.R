library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)


library(tibble)

df <- tribble(
  ~Survey, ~Predictor, ~Estimate, ~CI_low, ~CI_high, ~p,
  "Survey 1", "Male", -1.78, -3.97, 0.42, "0.112",
  "Survey 1", "18-29", -0.06, -3.58, 3.45, "0.972",
  "Survey 1", "30-44", -0.35, -3.32, 2.62, "0.815",
  "Survey 1", "45-59", 1.69, -1.93, 5.31, "0.357",
  "Survey 1", "≥60", -1.26, -4.40, 1.88, "0.429",
  
  "Survey 2", "Male", -6.67, -12.96, -0.38, "0.038",
  "Survey 2", "18-29", -3.44, -14.02, 7.14, "0.521",
  "Survey 2", "30-44", 16.95, 8.52, 25.39, "<0.001",
  "Survey 2", "45-59", 21.12, 8.45, 33.79, "0.001",
  "Survey 2", "≥60", 10.89, -0.61, 22.40, "0.063",
  "Survey 2", "SARS-CoV-2 IgG OD value in Survey  1", 5.56, 3.50, 7.62, "<0.001",
  "Survey 2", "1 Dose", 7.99, -1.03, 17.00, "0.082",
  
  "Survey 3", "Male", -5.80, -12.68, 1.09, "0.098",
  "Survey 3", "18-29", 10.72, -0.17, 21.61, "0.054",
  "Survey 3", "30-44", 3.40, -7.11, 13.91, "0.524",
  "Survey 3", "45-59", -6.40, -18.82, 6.01, "0.310",
  "Survey 3", "≥60", 11.88, -0.43, 24.19, "0.059",
  "Survey 3", "SARS-CoV-2 IgG OD value in Survey  1", 3.99, 1.48, 6.49, "0.002",
  "Survey 3", "SARS-CoV-2 IgG OD value in Survey  2", 1.81, 0.37, 3.25, "0.014",
  "Survey 3", "1 Dose", 5.44, -5.65, 16.52, "0.334",
  "Survey 3", "2 Doses", 21.50, 12.61, 30.38, "<0.001",
  "Survey 3", "3 Doses", 21.34, 11.90, 30.78, "<0.001",
  "Survey 3", "4 Doses", 28.01, 0.59, 55.44, "0.045",
  
  "Survey 4", "Male", -5.94, -12.34, 0.45, "0.068",
  "Survey 4", "18-29", -7.87, -17.76, 2.03, "0.119",
  "Survey 4", "30-44", -12.53, -22.25, -2.80, "0.012",
  "Survey 4", "45-59", -15.77, -26.65, -4.90, "0.005",
  "Survey 4", "≥60", -11.20, -22.64, 0.24, "0.055",
  "Survey 4", "SARS-CoV-2 IgG OD value in Survey  1", -0.04, -2.38, 2.31, "0.976",
  "Survey 4", "SARS-CoV-2 IgG OD value in Survey  2", 1.18, -0.16, 2.52, "0.084",
  "Survey 4", "SARS-CoV-2 IgG OD value in Survey  3", 2.84, 1.35, 4.33, "<0.001",
  "Survey 4", "1 Dose", -4.71, -16.26, 6.84, "0.422",
  "Survey 4", "2 Doses", -6.09, -15.82, 3.64, "0.219",
  "Survey 4", "3 Doses", 5.39, -4.07, 14.85, "0.262",
  "Survey 4", "4 Doses", -0.40, -10.03, 9.23, "0.935",
  
  "Survey 5", "Male", -1.01, -6.55, 4.54, "0.721",
  "Survey 5", "18-29", 1.36, -7.45, 10.17, "0.760",
  "Survey 5", "30-44", -10.78, -19.32, -2.25, "0.014",
  "Survey 5", "45-59", -8.26, -17.91, 1.40, "0.093",
  "Survey 5", "≥60", -1.50, -11.53, 8.52, "0.768",
  "Survey 5", "SARS-CoV-2 IgG OD value in Survey  1", 0.77, -1.28, 2.81, "0.459",
  "Survey 5", "SARS-CoV-2 IgG OD value in Survey  2", 0.63, -0.54, 1.80, "0.292",
  "Survey 5", "SARS-CoV-2 IgG OD value in Survey  3", 3.24, 1.92, 4.57, "<0.001",
  "Survey 5", "1 Dose", -6.32, -16.81, 4.17, "0.236",
  "Survey 5", "2 Doses", -4.81, -13.68, 4.06, "0.286",
  "Survey 5", "3 Doses", 2.08, -6.42, 10.57, "0.630",
  "Survey 5", "4 Doses", -3.94, -13.08, 5.21, "0.397",
  "Survey 5", "5 Doses", 8.80, -0.69, 18.28, "0.069"
)


library(dplyr)
library(ggplot2)
library(forcats)
library(patchwork)

# Calcular límites comunes del eje x
x_min <- floor(min(df$CI_low, na.rm = TRUE))
x_max <- ceiling(max(df$CI_high, na.rm = TRUE))

# Convertir p a numérico y crear columna de significancia
df <- df %>%
  mutate(
    p_num = suppressWarnings(as.numeric(p)),
    significant = ifelse((!is.na(p_num) & p_num < 0.05) | p == "<0.001", TRUE, FALSE)
  )


# Función para generar un forest plot por encuesta con orden invertido
plot_forest <- function(data, survey_title) {
  data <- data %>%
    mutate(Predictor = fct_rev(factor(Predictor, levels = unique(Predictor))))  # invertir orden
  
  ggplot(data, aes(x = Estimate, y = Predictor, color = significant)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "darkorange"), guide = "none") +
    xlim(x_min, x_max) +
    labs(
      title = paste("", survey_title),
      x = "OR (with 95% CI)",
      y = ""
    ) +
    theme_minimal(base_size = 14)
}

# Crear y mostrar todos los plots por encuesta
plots <- lapply(unique(df$Survey), function(s) {
  plot_forest(filter(df, Survey == s), s)
})

wrap_plots(plots, ncol = 3)




library(dplyr)
library(tidyr)

# Crear el data frame con columnas: survey, predictor, estimate, lower, upper, pvalue (todos como strings)
df <- tribble(
  ~Survey, ~Predictor, ~Estimate, ~CI_low, ~CI_high, ~p,
  
  "Survey 1", "Male", -0.04, -0.11, 0.03, "0.259",
  "Survey 1", "18-29", -0.04, -0.15, 0.07, "0.503",
  "Survey 1", "30-44", -0.01, -0.11, 0.08, "0.780",
  "Survey 1", "45-59", 0.03, -0.08, 0.15, "0.560",
  "Survey 1", "≥60", -0.07, -0.16, 0.03, "0.179",
  
  "Survey 2", "Male", -0.01, -0.10, 0.08, "0.827",
  "Survey 2", "18-29", -0.14, -0.29, 0.01, "0.071",
  "Survey 2", "30-44", 0.02, -0.10, 0.14, "0.773",
  "Survey 2", "45-59", 0.04, -0.14, 0.21, "0.689",
  "Survey 2", "≥60", 0.11, -0.05, 0.27, "0.182",
  "Survey 2", "SARS-CoV-2 IgG OD value in Survey 1", 0.04, 0.01, 0.07, "0.010",
  "Survey 2", "1 dose", 0.03, -0.09, 0.16, "0.603",

  "Survey 3", "Male", -0.01, -0.10, 0.08, "0.797",
  "Survey 3", "18-29", -0.01, -0.15, 0.13, "0.888",
  "Survey 3", "30-44", 0.03, -0.10, 0.16, "0.677",
  "Survey 3", "45-59", 0.16, -0.00, 0.31, "0.050",
  "Survey 3", "≥60", -0.07, -0.22, 0.09, "0.392",
  "Survey 3", "SARS-CoV-2 IgG OD value in Survey 1", -0.02, -0.05, 0.01, "0.167",
  "Survey 3", "SARS-CoV-2 IgG OD value in Survey 2", 0.02, 0.00, 0.04, "0.035",
  "Survey 3", "1 dose", -0.01, -0.15, 0.13, "0.902",
  "Survey 3", "2 doses", -0.01, -0.12, 0.10, "0.857",
  "Survey 3", "3 doses", -0.03, -0.15, 0.09, "0.662",
  "Survey 3", "4 doses", -0.12, -0.46, 0.23, "0.504",
  
  "Survey 4", "Male", 0.08, -0.02, 0.17, "0.103",
  "Survey 4", "18-29", -0.05, -0.20, 0.09, "0.479",
  "Survey 4", "30-44", -0.04, -0.19, 0.10, "0.554",
  "Survey 4", "45-59", 0.05, -0.11, 0.21, "0.518",
  "Survey 4", "≥60", -0.03, -0.20, 0.14, "0.750",
  "Survey 4", "SARS-CoV-2 IgG OD value in Survey 1", 0.01, -0.02, 0.05, "0.548",
  "Survey 4", "SARS-CoV-2 IgG OD value in Survey 2", 0.01, -0.01, 0.03, "0.376",
  "Survey 4", "SARS-CoV-2 IgG OD value in Survey 3", -0.02, -0.04, 0.01, "0.131",
  "Survey 4", "1 dose", 0.07, -0.10, 0.24, "0.404",
  "Survey 4", "2 doses", 0.08, -0.06, 0.23, "0.255",
  "Survey 4", "3 doses", -0.01, -0.15, 0.13, "0.874",
  "Survey 4", "4 doses", 0.07, -0.07, 0.21, "0.344",
 
  "Survey 5", "Male", 0.00, -0.08, 0.09, "0.930",
  "Survey 5", "18-29", -0.05, -0.19, 0.10, "0.526",
  "Survey 5", "30-44", -0.02, -0.16, 0.11, "0.722",
  "Survey 5", "45-59", -0.01, -0.16, 0.15, "0.912",
  "Survey 5", "≥60", -0.12, -0.28, 0.04, "0.152",
  "Survey 5", "SARS-CoV-2 IgG OD value in Survey 1", 0.00, -0.03, 0.03, "0.901",
  "Survey 5", "SARS-CoV-2 IgG OD value in Survey 2", 0.01, -0.01, 0.03, "0.167",
  "Survey 5", "SARS-CoV-2 IgG OD value in Survey 3", -0.03, -0.05, -0.01, "0.005",
  "Survey 5", "1 dose", 0.05, -0.12, 0.21, "0.580",
  "Survey 5", "2 doses", -0.01, -0.15, 0.13, "0.886",
  "Survey 5", "3 doses", 0.02, -0.11, 0.16, "0.718",
  "Survey 5", "4 doses", 0.06, -0.08, 0.21, "0.389",
  "Survey 5", "5 doses", -0.09, -0.24, 0.06, "0.225"
)



library(dplyr)
library(ggplot2)
library(forcats)
library(patchwork)

# Calcular límites comunes del eje x
x_min <- -0.5
x_max <- 0.5

# Convertir p a numérico y crear columna de significancia
df <- df %>%
  mutate(
    p_num = suppressWarnings(as.numeric(p)),
    significant = ifelse((!is.na(p_num) & p_num < 0.05) | p == "<0.001", TRUE, FALSE)
  )

# Función para generar un forest plot por encuesta con orden invertido
plot_forest <- function(data, survey_title) {
  data <- data %>%
    mutate(Predictor = fct_rev(factor(Predictor, levels = unique(Predictor))))  # invertir orden
  
  ggplot(data, aes(x = Estimate, y = Predictor, color = significant)) +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("FALSE" = "steelblue", "TRUE" = "darkorange"), guide = "none") +
    xlim(x_min, x_max) +
    labs(
      title = paste("", survey_title),
      x = "OR (with 95% CI)",
      y = ""
    ) +
    theme_minimal(base_size = 14)
}

# Crear y mostrar todos los plots por encuesta
plots <- lapply(unique(df$Survey), function(s) {
  plot_forest(filter(df, Survey == s), s)
})

wrap_plots(plots, ncol = 3)

