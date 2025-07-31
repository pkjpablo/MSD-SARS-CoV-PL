# Pablo Aguilar 
# IgG Re Vacunacion msd
# V1 

# paquetes 

library(readxl)
library(table1)
library(dplyr)
library(stringr)
library(tidyr)
library(gtsummary)

# Bancos ----
# Vacunacion 
db1 <- read_excel("dataset/bases anteriores/2024 06 05 Vacina covid resumen.xlsx")
names(db1)
# Sorologia
db2 <- read_excel("dataset/bases anteriores/2024 06 05 coleta elisa reinfec.xlsx")
names(db2)

#MSD
db5 <- read_excel("dataset/2025 02 25 205 participantsv2.xlsx")

# MSD 204 
db3 <- merge(db1,db2,by.x = "idnova", by.y = "id",all = T)

db5$pct.ancestral_p32.L49
db4 <- merge(db5[,c("idnova","pct.ancestral_p32.L49")],db3, by = "idnova", all.x=T)

db4 <- db4 %>% distinct(idnova, .keep_all = TRUE)

library(readr)
L49 <- read_csv("dataset/bases anteriores/L49Soroinquerito-PabloColetas_DATA_2025-01-17_0454.csv")
L49 <- L49[,c("idnova","dtcoleta")]

db4 <- merge(db4,L49, by = "idnova", all.x=T)
db4$dtcoletaL49 <- db4$dtcoleta

# limpiando  vacunacion ----
date_columns <- c("dtcoletaL45", "dtcoletaL46", "dtcoletaL47", "dtcoletaL48", "dtcoletaL49",
                  "data1", "data2", "data3", "data4", "data5", "data6")

db4[date_columns] <- lapply(db4[date_columns], as.Date)

count_effective_vaccines <- function(interview_date, vaccine_dates, effect_days = 14) {
  interview_date <- as.Date(interview_date) # Asegurarse de que interview_date sea Date
  valid_vaccines <- vaccine_dates[!is.na(vaccine_dates) & (vaccine_dates <= (interview_date - effect_days))]
  return(length(valid_vaccines))
}

# Aplicar la función para cada fecha de entrevista
db4$n.vacL45 <- apply(db4, 1, function(row) count_effective_vaccines(row["dtcoletaL45"], as.Date(c(row["data1"], row["data2"], row["data3"], row["data4"], row["data5"], row["data6"]))))
db4$n.vacL46 <- apply(db4, 1, function(row) count_effective_vaccines(row["dtcoletaL46"], as.Date(c(row["data1"], row["data2"], row["data3"], row["data4"], row["data5"], row["data6"]))))
db4$n.vacL47 <- apply(db4, 1, function(row) count_effective_vaccines(row["dtcoletaL47"], as.Date(c(row["data1"], row["data2"], row["data3"], row["data4"], row["data5"], row["data6"]))))
db4$n.vacL48 <- apply(db4, 1, function(row) count_effective_vaccines(row["dtcoletaL48"], as.Date(c(row["data1"], row["data2"], row["data3"], row["data4"], row["data5"], row["data6"]))))
db4$n.vacL49 <- apply(db4, 1, function(row) count_effective_vaccines(row["dtcoletaL49"], as.Date(c(row["data1"], row["data2"], row["data3"], row["data4"], row["data5"], row["data6"]))))

table(db4$n.vacL45)
table(db4$n.vacL46)
table(db4$n.vacL47)
table(db4$n.vacL48)
table(db4$n.vacL49)


db4$vac_cat.L45 <- ifelse(db4$n.vacL45 == 0, "no", "yes")
db4$vac_cat.L46<- ifelse(db4$n.vacL46 == 0, "no", "yes")
db4$vac_cat.L47 <- ifelse(db4$n.vacL47 == 0, "no", "yes")
db4$vac_cat.L48 <- ifelse(db4$n.vacL48 == 0, "no", "yes")
db4$vac_cat.L49 <- ifelse(db4$n.vacL49 == 0, "no", "yes")

# Limpiando previous exposure ----

db4$grupo <- paste(db4$razao_l45_grupo,
                   db4$razao_l46_grupo,
                   db4$RE,
                   db4$vac_cat.L46)

table1(~grupo, db4)

db4 <- db4 %>%
  mutate(grupo2 = recode(
    grupo,
    "NoSS1 NoSS2 NoRE no" = "NVSN",
    "NoSS1 NoSS2 NoRE yes" = "V",
    "NoSS1 SS2 NoRE no" = "SS2",
    "NoSS1 SS2 NoRE yes" = "SS2 or V",
    "SS1 NoSS2 NoRE no" = "SS1",
    "SS1 NoSS2 NoRE V" = "SV",
    "SS1 SS2 NoRE no" = "SS1",
    "SS1 SS2 NoRE yes" = "SV",
    "SS1 SS2 RE no" = "RE",
    "SS1 SS2 RE yes" = "RE or SV"
  ))

table1(~grupo2, db4)

db4$grupo2 <- factor(db4$grupo2, levels = c("NVSN",
                                            "V",
                                            "SS1",
                                            "SS2",
                                            "RE",
                                            "SV",
                                            "SS2 or V",
                                            "RE or SV"))
table1(~grupo2, db4)


#XXXXXXX
L45socio <- read_excel("dataset/bases anteriores/L45 sociodemografico pablo.xlsx")

socio <- L45socio %>%
  select(idnova,cens_idade, idade_cat, cind_sexo, raca2, serie, schooling, estado_civil3, trabalho2)

dbx <- merge(socio, db4,by = "idnova",all.y = T)
dbx$categoria_edad <- cut(dbx$cens_idade,
                          breaks = c(-Inf, 5, 18, 34, 64, Inf),
                          labels = c("≤5", "6-18", "19-34", "34-64", ">64"),
                          right = TRUE)

dbx$cens_idade
dbx$razao_l45_grupo
#dbx$razao_grupo.L45
tab_socio <- dbx%>% 
  select(cens_idade,categoria_edad,cind_sexo,
         vac_cat.L45,vac_cat.L46,vac_cat.L47,vac_cat.L48,vac_cat.L49,
         razao_l45_grupo,
         razao_l46_grupo,
         razao_l47_grupo,
         grupo2) %>% 
  tbl_summary()

tab_socio

dbx$new <- ifelse(is.na(dbx$pct.ancestral_p32.L49), "No", "Yes")

tab_socio <- dbx%>% 
  select(cens_idade,categoria_edad,cind_sexo,
         vac_cat.L45,vac_cat.L46,vac_cat.L47,vac_cat.L48,vac_cat.L49,
         razao_l45_grupo,
         razao_l46_grupo,
         razao_l47_grupo,
         grupo2,new) %>% 
  tbl_summary(by = new)
tab_socio

#library(writexl)
write_xlsx(dbx, 'dataset/bases anteriores/socio2.xlsx')


library(readxl)
dbx <- read_excel("dataset/bases anteriores/socio2.xlsx")
dbx$razao_l45_grupo

dbx$grupo3 <- case_when(
  dbx$grupo2 %in% c("NVSN") ~ "0expo",
  dbx$grupo2 %in% c("V", "SS1", "SS2", "SS2 or V") ~ "1expo",
  dbx$grupo2 %in% c("RE", "SV", "RE or SV") ~ "2expo",
  TRUE ~ NA_character_  # Asigna NA a cualquier valor que no esté en las categorías definidas
)

tab_socio <- dbx%>% 
  select(cens_idade,categoria_edad,cind_sexo,
         vac_cat.L45,vac_cat.L46,vac_cat.L47,vac_cat.L48,
         razao_l45_grupo,
         razao_l46_grupo,
         razao_l47_grupo,
         grupo3) %>% 
  tbl_summary()

tab_socio



# Transformar de formato ancho a largo
db4_long <- db4 %>%
  pivot_longer(
    cols = c(starts_with("dtcoleta"),
             starts_with("razao_"),
             starts_with("razao_grupo."),
             starts_with("vacL"),
             starts_with("vac_cat."),
             starts_with("n.vac")),
    names_to = c(".value", "survey"),
    names_pattern = "(.*)(L\\d+)$"
  )


#write_xlsx(db4_long, 'dataset/bases anteriores/2024 06 05 137 igg vac long.xlsx')

library(readxl)
L45socio <- read_excel("dataset/bases anteriores/L45 sociodemografico pablo.xlsx")

socio <- L45socio %>%
  select(idnova,cens_idade, idade_cat, cind_sexo, raca2, serie, schooling, estado_civil3, trabalho2)

db6 <- merge(socio, db4_long,by = "idnova",all.y = T)


#write_xlsx(db4_long, 'dataset/bases anteriores/2024 06 05 137 igg vac + socio long.xlsx')

## MSD
#msd <- read_excel("dataset/2024 06 05 205 participantsv2.xlsx")
msd <- read_excel("dataset/2025 02 25 205 participantsv2.xlsx")

#msd$`pct.ancestral.L49_Panel 32` <- msd$`pct.ancestral.L49_Panel 33`

names(msd)

msd_long <- msd %>%
  pivot_longer(
    cols = starts_with("pct."),
    names_to = c("lineage", "survey"),
    names_pattern = "^pct\\.(.*)\\.L(\\d+)$",
    values_to = "value"
  ) %>%
  mutate(
    panel = case_when(
      lineage %in% c("ancestral_p23", "AY.4.2", "B.1.1.7", "B.1.351", "B.1.617.Seq1", "B.1.617.Seq2",  "P.1","BA.1_p23") ~ 23,
      lineage %in% c("ancestral_p32", "BA.1_p32", "BA.2.75", "BA.2.75.2", "BA.4.6", "BA.5", "BF.7", "BQ.1", "BQ.1.1", "XBB.1") ~ 32,
      TRUE ~ NA_real_  # Si no está en las listas, se asigna NA
    )
  )


msd_long$lineage <- as.factor(msd_long$lineage)
levels(msd_long$lineage)

msd_long <- subset(msd_long, lineage %in% c(
  "ancestral_p23", "B.1.1.7", "B.1.351", "B.1.617.Seq2", 
  "BA.1_p32", "BA.2.75", "BA.5", "BF.7", "BQ.1", 
  "BQ.1.1", "P.1", "XBB.1"
))

msd_long$lineage <- gsub("_p23", "", msd_long$lineage)  # Elimina "_p23"
msd_long$lineage <- gsub("_p32", "", msd_long$lineage)  # Elimina "_p32"

# Actualizar los niveles de la columna 'lineage'
msd_long$lineage <- factor(msd_long$lineage)

# Verificar los niveles actualizados
levels(msd_long$lineage)

# todo junto 

db6$survey
db6$idnova

msd_long$survey
msd_long$idnova

msd_long <- msd_long %>%
  mutate(
    survey = if_else(
      !is.na(survey),  # Si survey no es NA
      paste0("L", survey),  # Agrega "L" al valor de survey
      survey  # Si es NA, no cambia
    )
  )

db6$survey
msd_long$survey

db7 <- merge(db6, msd_long, by = c("idnova", "survey"))
db7$lineage <- as.factor(db7$lineage)
# Define el nuevo orden de las variantes basado en su aparición
nuevo_orden <- c("ancestral", "B.1.1.7", "B.1.351", "P.1", 
                 "B.1.617.Seq1", "B.1.617.Seq2", "AY.4.2", "BA.1", 
                 "BA.2.75", "BA.2.75.2", "BA.5", "BA.4.6", 
                 "BF.7", "BQ.1", "BQ.1.1", "XBB.1")

# Reordena los niveles de la variable 'lineage'
db7$lineage <- factor(db7$lineage, levels = nuevo_orden)

db7$n.vac

library(ggplot2)
hex <- c("#FF0000", "#FFA500", "#FFFF00", "#008000", "#9999ff", "#000066")

# Asumiendo que 'db7' es tu dataframe y que 'panel', 'survey', 'value', 'idnova', 'grupo2', 'lineage' y 'n.vac' son columnas en 'db7'


db7$grupo2 <- factor(db7$grupo2, levels = c("NVSN",
                                            "V",
                                            "SS2",
                                            "SS2 or V",
                                            "SS1",
                                            "RE",
                                            "SV",
                                            "RE or SV"))
table1(~grupo2, db7)

db7$grupo3 <- case_when(
  db7$grupo2 %in% c("NVSN") ~ "0expo",
  db7$grupo2 %in% c("V", "SS1", "SS2", "SS2 or V") ~ "1expo",
  db7$grupo2 %in% c("RE", "SV", "RE or SV") ~ "2expo",
  TRUE ~ NA_character_  # Asigna NA a cualquier valor que no esté en las categorías definidas
)

# Convertir a factor con niveles ordenados
db7$grupo3 <- factor(db7$grupo3, levels = c("0expo", "1expo", "2expo"))

table1(~grupo2, db7)

db7$grupo2 <- factor(db7$grupo2, levels = c("NVSN", "V", "SS2", "SS2 or V", "SS1", "RE", "SV", "RE or SV"))
db7$lineage <- factor(db7$lineage)

ggplot(subset(db7, panel == "23"), aes(x = survey, y = value, group = idnova, color = n.vac)) +
  geom_line() +
  geom_point() +  # Añade puntos al gráfico
  scale_color_gradientn(colours = rev(hex)) +  # Ajusta el número dentro de `hex` según sea necesario
  facet_grid(grupo2 ~ lineage) +
  theme_minimal() +
  labs(x = "Survey", y = "% inhibition", title = "Panel 23") +
  theme(legend.position = "right")

ggplot(subset(db7, panel == "32"), aes(x = survey, y = value, group = idnova, color = n.vac)) +
  geom_line() +
  geom_point() +  # Añade puntos al gráfico
  scale_color_gradientn(colours = rev(hex)) +  # Ajusta el número dentro de `hex` según sea necesario
  facet_grid(grupo2 ~ lineage) +
  theme_minimal() +
  labs(x = "Survey", y = "% inhibition", title = "Panel 32") +
  theme(legend.position = "right")


## 3 grupos 
ggplot(subset(db7, panel == "23"), aes(x = survey, y = value, group = idnova, color = n.vac)) +
  geom_line() +
  geom_point() +  # Añade puntos al gráfico
  scale_color_gradientn(colours = rev(hex)) +  # Ajusta el número dentro de `hex` según sea necesario
  facet_grid(grupo3 ~ lineage) +
  theme_minimal() +
  labs(x = "Survey", y = "% inhibition", title = "Panel 23") +
  theme(legend.position = "right")

ggplot(subset(db7, panel == "32"), aes(x = survey, y = value, group = idnova, color = n.vac)) +
  geom_line() +
  geom_point() +  # Añade puntos al gráfico
  scale_color_gradientn(colours = rev(hex)) +  # Ajusta el número dentro de `hex` según sea necesario
  facet_grid(grupo3 ~ lineage) +
  theme_minimal() +
  labs(x = "Survey", y = "% inhibition", title = "Panel 32") +
  theme(legend.position = "right")


## boxplot 
medianas <- subset(db7, panel == "23") %>%
  group_by(grupo3, lineage, survey) %>%
  summarize(mediana = median(value, na.rm = TRUE), .groups = "drop")

ggplot(subset(db7, panel == "23"), aes(x = survey, y = value, color = n.vac)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(fill = survey)) +
  geom_line(
    data = medianas,
    aes(x = survey, y = mediana, group = interaction(grupo3, lineage)),
    color = "black", size = 1
  ) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.2) +
  scale_color_gradientn(colours = rev(hex)) +
  scale_fill_manual(values = c("#FFDDC1", "#FFABAB", "#FFC3A0", "#D5AAFF","#FFC3A9")) +
  facet_grid(grupo3 ~ lineage) +
  theme_minimal() +
  labs(
    x = "Survey",
    y = "% Inhibition",
    title = "Panel 23"
  ) 


medianas <- subset(db7, panel == "32") %>%
  group_by(grupo3, lineage, survey) %>%
  summarize(mediana = median(value, na.rm = TRUE), .groups = "drop")

ggplot(subset(db7, panel == "32"), aes(x = survey, y = value, color = n.vac)) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, aes(fill = survey)) +
  geom_line(
    data = medianas,
    aes(x = survey, y = mediana, group = interaction(grupo3, lineage)),
    color = "black", size = 1
  ) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.2) +
  scale_color_gradientn(colours = rev(hex)) +
  scale_fill_manual(values = c("#FFDDC1", "#FFABAB", "#FFC3A0", "#D5AAFF","#FFAAFF")) +
  facet_grid(grupo3 ~ lineage) +
  theme_minimal() +
  labs(
    x = "Survey",
    y = "% Inhibition",
    title = "Panel 23"
  ) 




library(dplyr)
library(rstatix)

# Filtrar el subconjunto relevante
subset_db7 <- subset(db7, panel == "23")

# Realizar las comparaciones estadísticas
signif_table <- subset_db7 %>%
  group_by(grupo3, lineage) %>%  # Agrupar por grupo3 y lineage
  pairwise_t_test(
    value ~ survey,  # Comparación de `value` entre niveles de `survey`
    p.adjust.method = "bonferroni"  # Ajuste para comparaciones múltiples
  ) %>%
  ungroup()

# Mostrar tabla de significancia
signif_table


library(dplyr)
library(rstatix)

# Filtrar el subconjunto relevante
subset_db7 <- subset(db7, panel == "23")

# Comparar niveles consecutivos de `survey`
simplified_signif_table <- subset_db7 %>%
  group_by(grupo3, lineage) %>%  # Agrupar por grupo3 y lineage
  filter(!is.na(survey)) %>%
  arrange(survey) %>%
  mutate(survey = factor(survey, levels = unique(survey))) %>%
  pairwise_t_test(
    value ~ survey,  # Comparar `value` entre niveles de `survey`
    p.adjust.method = "none",  # Sin ajuste de p-value
    comparisons = list(
      c("L45", "L46"),
      c("L46", "L47"),
      c("L47", "L48")
    )  # Comparaciones específicas
  ) %>%
  ungroup()

# Mostrar tabla simplificada
simplified_signif_table


# Filtrar el subconjunto relevante
subset_db7 <- subset(db7, panel == "32")

# Comparar niveles consecutivos de `survey`
simplified_signif_table <- subset_db7 %>%
  group_by(grupo3, lineage) %>%  # Agrupar por grupo3 y lineage
  filter(!is.na(survey)) %>%
  arrange(survey) %>%
  mutate(survey = factor(survey, levels = unique(survey))) %>%
  pairwise_t_test(
    value ~ survey,  # Comparar `value` entre niveles de `survey`
    p.adjust.method = "none",  # Sin ajuste de p-value
    comparisons = list(
      c("L45", "L46"),
      c("L46", "L47"),
      c("L47", "L48"),c("L48", "L49")
    )  # Comparaciones específicas
  ) %>%
  ungroup()

# Mostrar tabla simplificada
simplified_signif_table





library(writexl)
write_xlsx(db7, 'dataset/bases anteriores/2024 06 05 192 long primiera parte.xlsx')

#write_xlsx(db7, 'dataset/bases anteriores/2024 06 05 204 long primiera parte.xlsx')

#write_xlsx(db7, 'dataset/bases anteriores/2024 06 05 137 long primiera parte.xlsx')


db7 <- db7 %>%
  mutate(grupo2 = case_when(
    survey == "L45" & grupo2 %in% c("NVSN", "V", "SS2", "SS2 or V") ~ "Negative",
    survey == "L45" & grupo2 %in% c("SS1", "SV", "RE", "RE or SV") ~ "Positive",
    TRUE ~ grupo2
  ))
db7

db7 <- db7 %>%
  mutate(variant = paste(lineage, panel, sep = "_"))

order_variants <- c(
  "ancestral_23", "ancestral_32",
  "B.1.1.7_23",
  "B.1.351_23",
  "P.1_23",
  "B.1.617.Seq1_23",
  "B.1.617.Seq2_23",
  "AY.4.2_23", 
  "BA.1_23", "BA.1_32",
  "BA.2.75_32",
  "BA.2.75.2_32",
  "BA.4.6_32",
  "BA.5_32",
  "BF.7_32",
  "BQ.1.1_32", "BQ.1_32",
  "XBB.1_32"
)

# Reordenar la variable en el dataframe
db7$variant <- factor(db7$variant, levels = order_variants)

table(db7$variant)

db7$n.vac

## Panel 23

# L45 

db7$grupo2 <- factor(db7$grupo2, levels = c("NVSN",
                                            "V",
                                            "SS1",
                                            "SS2",
                                            "SS2 or V",
                                            "SV",
                                            "RE",
                                            "RE or SV"))


a2 <- ggplot(subset(db7, survey == "L45"), aes(x = variant, y = value, color = lineage)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +  
  stat_summary(fun = mean, geom = "point", shape = 1, fill = "black", size = 3) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid", size = 1) +
  facet_grid(rows = vars(grupo2)) +  # Eliminar categorías vacías
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

a2



b2 <- ggplot(subset(db7, survey == "L46"), aes(x = variant, y = value, color = lineage)) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +  
  stat_summary(fun = mean, geom = "point", shape = 1, fill = "black", size = 3) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid", size = 1) +
  facet_grid(rows = vars(grupo2)) +  # Eliminar categorías vacías
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

b2


c2 <- ggplot(subset(db7, survey == "L47"), aes(x = variant, y = value, color = factor(n.vac))) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +  
  ylim(-0.5, 1.2) +
  stat_summary(fun = mean, geom = "point", shape = 1, fill = "black", size = 0.1) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid", size = 0.5) +
  facet_grid(rows = vars(grupo2)) +
  scale_color_discrete(name = "Dosis de Vacuna") +  # Opcional: para personalizar la leyenda
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

c2


d2 <- ggplot(subset(db7, survey == "L48"), aes(x = variant, y = value, color = factor(n.vac))) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +  
  ylim(-0.5, 1.2) +
  stat_summary(fun = mean, geom = "point", shape = 1, fill = "black", size = 0.1) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid", size = 0.5) +
  facet_grid(rows = vars(grupo2)) +
  scale_color_discrete(name = "Dosis de Vacuna") +  # Opcional: para personalizar la leyenda
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

d2

d2 <- ggplot(subset(db7, survey == "L49"), aes(x = variant, y = value, color = factor(n.vac))) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.5) +  
  ylim(-0.5, 1.2) +
  stat_summary(fun = mean, geom = "point", shape = 1, fill = "black", size = 0.1) +
  stat_summary(fun.data = mean_sdl, fun.args = list(mult = 1), geom = "errorbar", width = 0.2, color = "black") +
  stat_summary(fun = mean, geom = "line", aes(group = 1), color = "black", size = 1) +
  geom_hline(yintercept = 0, color = "darkred", linetype = "solid", size = 0.5) +
  facet_grid(rows = vars(grupo2)) +
  scale_color_discrete(name = "Dosis de Vacuna") +  # Opcional: para personalizar la leyenda
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

d2





