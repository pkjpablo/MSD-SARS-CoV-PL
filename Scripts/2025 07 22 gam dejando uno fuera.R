# Pablo 
# SARS-CoV-2 Landscape GAM
# surveys 


## Libraries ----

library(dplyr)
library(mgcv)
library(plotly)
library(stringr)



## Dataset----
#saveRDS(df_coords, file = "dataset/df_coords.rds")
df_coords <- readRDS("dataset/df_coords.rds")

#saveRDS(db_long, file = "dataset/db_long.rds")
db_long <- readRDS("dataset/db_long.rds")

db_long ## MSD Pau da Lima 
df_coords ## Roesler cordenadas 


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


table(db_long$ag_name)
table(df_coords$ag_name)
mapColors


evaluate_variant_prediction_v2 <- function(survey_name, db_long, df_coords) {
  # Subset survey
  db_sub <- db_long %>% filter(str_starts(sr_name, survey_name))
  merged <- merge(db_sub, df_coords, by = "ag_name", all.x = TRUE)
  
  # Obtener promedio z por variante y posición (x,y)
  data_avg <- merged %>%
    mutate(z = as.numeric(titer)) %>%
    group_by(variant = ag_name, x, y) %>%
    summarise(z = mean(z, na.rm = TRUE), .groups = "drop")
  
  # GAM completo con todas las variantes (promedio por variante)
  gam_full <- gam(z ~ s(x, y, k = 11), data = data_avg)
  
  results <- list()
  
  variants <- unique(data_avg$variant)
  
  for (v in variants) {
    # posición (x,y) de la variante
    pos_var <- data_avg %>% filter(variant == v) %>% select(x, y, z)
    
    # Predicción con GAM completo (modelo con todos)
    pred_full <- predict(gam_full, newdata = pos_var)
    
    # Excluir la variante completa del dataset para ajustar nuevo GAM
    data_train <- data_avg %>% filter(variant != v)
    gam_excl <- gam(z ~ s(x, y, k = 11), data = data_train)
    
    # Predicción en la posición de la variante excluida
    pred_excl <- predict(gam_excl, newdata = pos_var)
    
    df_res <- data.frame(
      Survey = survey_name,
      Variant = v,
      Real_z = pos_var$z,
      Predicted_z_full = pred_full,
      Predicted_z_excl = pred_excl
    )
    
    df_res$Abs_Error <- abs(df_res$Predicted_z_excl - df_res$Real_z)
    df_res$Rel_Error <- df_res$Abs_Error / df_res$Real_z * 100
    
    results[[v]] <- df_res
  }
  
  do.call(rbind, results)
}



result_L45 <- evaluate_variant_prediction_v2("L45", db_long, df_coords)
head(result_L45)


surveys <- c("L45", "L46", "L47", "L48", "L49")


all_results <- lapply(surveys, function(sv) evaluate_variant_prediction_v2(sv, db_long, df_coords))
combined_results <- do.call(rbind, all_results)

print(combined_results)

library(dplyr)

# Define el orden deseado de las variantes
variant_order <- c(
  "Ancestral", "Alpha", "Beta", "Delta", "Gamma",
  "BA.1", "BA.2", "BA.5", "BF.7", "BQ.1", "BQ.1.1",
  "XBB.1"
)

# Combina la lista en un único data.frame
combined_results <- bind_rows(all_results)

# Ordena según el orden deseado
combined_results <- combined_results %>%
  mutate(Variant = factor(Variant, levels = variant_order)) %>%
  arrange(Variant)

# Imprimir resultado
print(combined_results)


library(ggplot2)

ggplot(combined_results, aes(x = reorder(Variant, -Rel_Error), y = Rel_Error)) +
  geom_boxplot() +
  coord_flip() +
  labs(title = "Errores relativos por variante",
       y = "Error relativo (%)", x = "Variante") +
  theme_minimal()

library(tidyr)
library(dplyr)

combined_long <- combined_results %>%
  pivot_longer(cols = c(Real_z, Predicted_z_full, Predicted_z_excl),
               names_to = "Tipo",
               values_to = "Valor")

library(ggplot2)

ggplot(combined_long, aes(x = Variant, y = Valor, color = Tipo, shape = Tipo)) +
  geom_point(size = 3, position = position_dodge(width = 0.5)) +
  facet_wrap(~Survey, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Valor real vs predicciones por variante y survey",
       x = "Variante", y = "Valor Z normalizado")

