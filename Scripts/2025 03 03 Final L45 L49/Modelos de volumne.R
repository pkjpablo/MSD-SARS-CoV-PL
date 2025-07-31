# Modelos landscape
# Pablo

## use the version 4.4.1 

#install.packages("sjPlot")
# libraries
library(readxl)
library(lme4)
library(sjPlot)
library(performance)
library(splines)
library(forcats)


setwd("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/1. Landscape clean/MSD")

df45 <- readRDS(file = "dataset/volumenes/L45_volumenx2.rds")
df46 <- readRDS(file = "dataset/volumenes/L46_volumenx2.rds")
df47 <- readRDS(file = "dataset/volumenes/L47_volumenx2.rds")
df48 <- readRDS(file = "dataset/volumenes/L48_volumenx2.rds")
df49 <- readRDS(file = "dataset/volumenes/L49_volumenx2.rds")

file_path <- "dataset/volumenes/L45_49_all2.rds"

if (!file.exists(file_path)) {
  db01 <- rbind(df45, df46, df47, df48, df49)
  saveRDS(db01, file = file_path)
} else {
  db01 <- readRDS(file_path)
}


L45=read.csv('dataset/bases anteriores/L45SoroinqueritoHarm-L45Basico_DATA_2024-03-07_1402.csv')
L45 <- L45[, c("idnova","cind_sexo", "idade_calc_check", "razao_l45", "dtcoleta")]
cols_to_modify <- setdiff(names(L45), c("idnova", "razao_l45"))
names(L45)[names(L45) %in% cols_to_modify] <- paste0(cols_to_modify, "_L45")

L46=read.csv('dataset/bases anteriores/L46IFAASHarmonizado-L46Basico_DATA_2024-03-07_1403.csv')
L46 <- L46[, c("idnova","cind_sexo", "idade_calc_check", "razao_l46", "dtcoleta")]
cols_to_modify <- setdiff(names(L46), c("idnova", "razao_l46"))
names(L46)[names(L46) %in% cols_to_modify] <- paste0(cols_to_modify, "_L46")

L47=read.csv('dataset/bases anteriores/L47SoroinqueritoHarm-PabloBasico_DATA_2025-03-06_1413.csv')
L47 <- L47[, c("idnova","cind_sexo", "idade_calc_check", "razao_l47", "dtcoleta")]
cols_to_modify <- setdiff(names(L47), c("idnova", "razao_l47"))
names(L47)[names(L47) %in% cols_to_modify] <- paste0(cols_to_modify, "_L47")

L48=read.csv('dataset/bases anteriores/L48Soroinquerito-PabloBasico_DATA_2025-03-06_1404.csv')
L48$razao_l48 <- NA
L48 <- L48[, c("idnova","cind_sexo", "idade_calc_check", "razao_l48", "dtcoleta")]
cols_to_modify <- setdiff(names(L48), c("idnova", "razao_l48"))
names(L48)[names(L48) %in% cols_to_modify] <- paste0(cols_to_modify, "_L48")

L49=read.csv('dataset/bases anteriores/L49Soroinquerito-PabloBasico_DATA_2025-03-14_1113.csv')
L49$razao_l49 <- NA
L49 <- L49[, c("idnova","cind_sexo", "idade_calc_check", "razao_l49", "dtcoleta")]
cols_to_modify <- setdiff(names(L49), c("idnova", "razao_l49"))
names(L49)[names(L49) %in% cols_to_modify] <- paste0(cols_to_modify, "_L49")


# Lista de dataframes con nombres
dfs <- list(L45 = L45, L46 = L46, L47 = L47, L48 = L48, L49 = L49)

dfs <- lapply(dfs, function(df) df[!is.na(df$idnova), ])

# Función para fusionar con sufijos adecuados
merge_with_suffix <- function(x, y, name_y) {
  merge(x, y, by = "idnova", all = TRUE)
}

# Aplicar merge secuencialmente con nombres correctos
db_merged <- Reduce(function(x, y) {
  merge_with_suffix(x, y, deparse(substitute(y)))
}, dfs)


db01$idnova <- db01$id
db01$idnova <- substr(db01$idnova, 5, nchar(db01$idnova))

idnova_unicos <- unique(db01$idnova)
db_merged_filtrado <- db_merged[db_merged$idnova %in% idnova_unicos, ]
db_merged_filtrado <- subset(db_merged_filtrado, !is.na(razao_l46))

# Identificar los idnova duplicados
duplicados_idnova <- db_merged_filtrado$idnova[duplicated(db_merged_filtrado$idnova) | duplicated(db_merged_filtrado$idnova, fromLast = TRUE)]

# Crear un nuevo dataset con las filas duplicadas
db_duplicados <- db_merged_filtrado[db_merged_filtrado$idnova %in% duplicados_idnova, ]



library(dplyr)

db_merged_filtrado <- db_merged_filtrado %>%
  distinct(idnova, .keep_all = TRUE)



names(db_merged_filtrado)

names(db_merged_filtrado) <- c(
  "idnova", "cind_sexo_l45", "idade_calc_check_l45", "razao_l45", "dtcoleta_l45",
  "cind_sexo_l46", "idade_calc_check_l46", "razao_l46", "dtcoleta_l46",
  "cind_sexo_l47", "idade_calc_check_l47", "razao_l47", "dtcoleta_l47",
  "cind_sexo_l48", "idade_calc_check_l48", "razao_l48", "dtcoleta_l48",
  "cind_sexo_l49", "idade_calc_check_l49", "razao_l49", "dtcoleta_l49"
)

summary(db_merged_filtrado)


names(db_merged_filtrado)






### limpiando vacunacion 

db1 <- read_excel("dataset/bases anteriores/2024 06 05 Vacina covid resumen.xlsx")

db4 <- merge(db_merged_filtrado,db1,by = "idnova",all.x = T)
names(db4)

date_columns <- c("dtcoleta_l45", "dtcoleta_l46", "dtcoleta_l47", "dtcoleta_l48", "dtcoleta_l49",
                  "data1", "data2", "data3", "data4", "data5", "data6")

db4[date_columns] <- lapply(db4[date_columns], as.Date)

count_effective_vaccines <- function(interview_date, vaccine_dates, effect_days = 14) {
  interview_date <- as.Date(interview_date) # Asegurarse de que interview_date sea Date
  valid_vaccines <- vaccine_dates[!is.na(vaccine_dates) & (vaccine_dates <= (interview_date - effect_days))]
  return(length(valid_vaccines))
}

# Aplicar la función para cada fecha de entrevista
db4$n.vacL45 <- apply(db4, 1, function(row) count_effective_vaccines(row["dtcoleta_l45"], as.Date(c(row["data1"], row["data2"], row["data3"], row["data4"], row["data5"], row["data6"]))))
db4$n.vacL46 <- apply(db4, 1, function(row) count_effective_vaccines(row["dtcoleta_l46"], as.Date(c(row["data1"], row["data2"], row["data3"], row["data4"], row["data5"], row["data6"]))))
db4$n.vacL47 <- apply(db4, 1, function(row) count_effective_vaccines(row["dtcoleta_l47"], as.Date(c(row["data1"], row["data2"], row["data3"], row["data4"], row["data5"], row["data6"]))))
db4$n.vacL48 <- apply(db4, 1, function(row) count_effective_vaccines(row["dtcoleta_l48"], as.Date(c(row["data1"], row["data2"], row["data3"], row["data4"], row["data5"], row["data6"]))))
db4$n.vacL49 <- apply(db4, 1, function(row) count_effective_vaccines(row["dtcoleta_l49"], as.Date(c(row["data1"], row["data2"], row["data3"], row["data4"], row["data5"], row["data6"]))))

18+35+7
table(db4$n.vacL45,useNA = "always")
table(db4$n.vacL46,useNA = "always")
table(db4$n.vacL47,useNA = "always")
table(db4$n.vacL48,useNA = "always")
table(db4$n.vacL49,useNA = "always")


db4$vac_cat.L45 <- ifelse(db4$n.vacL45 == 0, "no", "yes")
db4$vac_cat.L46<- ifelse(db4$n.vacL46 == 0, "no", "yes")
db4$vac_cat.L47 <- ifelse(db4$n.vacL47 == 0, "no", "yes")
db4$vac_cat.L48 <- ifelse(db4$n.vacL48 == 0, "no", "yes")
db4$vac_cat.L49 <- ifelse(db4$n.vacL49 == 0, "no", "yes")


#### Elisa 

db2 <- read_excel("dataset/bases anteriores/2024 06 05 coleta elisa reinfec.xlsx")
names(db2)
names(db4)

db2_selected <- db2[, c("id", "razao_l45_grupo", "razao_l46_grupo", "RE")]


db5 <- merge(db4,db2_selected,by.x = "idnova",by.y = "id",all.x = T)

duplicados_idnova <- db5$idnova[duplicated(db5$idnova) | duplicated(db5$idnova, fromLast = TRUE)]
# Crear un nuevo subset con las filas que contienen los idnova duplicados
db5_duplicados <- db5[db5$idnova %in% duplicados_idnova, ]

db4 <- unique(db5)

db4$grupo <- paste(db4$razao_l45_grupo,
                   db4$razao_l46_grupo,
                   db4$RE,
                   db4$vac_cat.L46)

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

table1::table1(~grupo2, db4)

db4$grupo2 <- factor(db4$grupo2, levels = c("NVSN",
                                            "V",
                                            "SS1",
                                            "SS2",
                                            "RE",
                                            "SV",
                                            "SS2 or V",
                                            "RE or SV"))

db4$grupo3L46 <- case_when(
  db4$grupo2 %in% c("NVSN") ~ "0expo",
  db4$grupo2 %in% c("V", "SS1", "SS2", "SS2 or V") ~ "1expo",
  db4$grupo2 %in% c("RE", "SV", "RE or SV") ~ "2expo",
  TRUE ~ NA_character_  # Asigna NA a cualquier valor que no esté en las categorías definidas
)


names(db4)
names(db4) <- gsub("(l45|l46|l47|l48|l49)", "\\U\\1", names(db4), perl = TRUE)

db4$n.vacL45

# Cargar librería tidyr
library(tidyr)

# Seleccionar las columnas relacionadas con los surveys
survey_columns <- grep("L45|L46|L47|L48|L49", names(db4), value = TRUE)


names(db4)
db4_long <- db4 %>%
  pivot_longer(
    cols = matches("L(45|46|47|48|49)$"),  # columnas que terminan en L45, L46, ...
    names_to = c(".value", "survey"),     # divide en valor base + número del lote
    names_pattern = "(.*)_L(\\d+)"        # patrón corregido con L mayúscula
  ) %>%
  mutate(idnova_survey = paste(idnova, survey, sep = "_"))

db4 <- db4[, !(names(db4) %in% c("razao_L45_grupo", "razao_L46_grupo"))]
names(db4)

# Convertir el dataframe de wide a long
db4_long <- db4 %>%
  pivot_longer(
    cols = matches("L(45|46|47|48|49)"),  # Selecciona todas las columnas con L45-L49 en cualquier parte
    names_to = c(".value", "survey"),  
    names_pattern = "(.*)_L?(\\d+)"  # Extrae la parte base y el número de la encuesta, manejando "L" opcionalmente
  ) %>%
  mutate(idnova_survey = paste(idnova, survey, sep = "_"))  # Crear nueva ID única


# Verificar los primeros registros del dataframe long
head(db4_long)


db4_long$idjp <- paste0(db4_long$survey,"_",db4_long$idnova)
db4_long$idjp <- paste0("L", db4_long$idjp)
db4_long$survey <- paste0("L", db4_long$survey)

db01$volumen100 <- db01$volumen/270*100



db_v <- merge(db01,db4_long, by.x = "id",by.y = "idjp",all.x = T)

db_v$categoria_edad <- cut(db_v$idade_calc_check,
                            breaks = c(-Inf, 18, 30, 45, 59, Inf),
                            labels = c("<18", "18-29", "30-44", "45-59", "≥60"),
                            right = TRUE)



db_v$cind_sexo <- factor(db_v$cind_sexo, levels = c(0, 1), labels = c("female", "male"))


db_v <- merge(db_v,db4,by.x ="idnova.x",by.y ="idnova",all.x = T)


db45 <- subset(db_v, survey=="L45")
db46 <- subset(db_v, survey=="L46")
db47 <- subset(db_v, survey=="L47")
db48 <- subset(db_v, survey=="L48")
db49 <- subset(db_v, survey=="L49")

db45$ce

summary(db45)

############ Modelo L45 

db45 <- db45 %>%
  left_join(
    dfL45.2 %>% select(id, volumen_percent, slope_gam),
    by = c("id" = "id")
  )



db45$volumen100
db45$volumen10

m45 <- lm(volumen100~ 
            cind_sexo +
            categoria_edad+razao_l45_grupo  , data = db45)

tab_model(m45,show.aic = T)

backwards = step(m45)

tab_model(backwards,show.aic = T)

m45 <- lm(volumen100~ 
            cind_sexo +
            categoria_edad  , data = db45)

tab_model(m45,show.aic = T)


############ Modelo L46
db46$vac_cat.L46.x
m46 <- lm(volumen100 ~  cind_sexo+
            categoria_edad+razao_l45+factor(vac_cat.L46.x), data = db46)

tab_model(m46)

library(dplyr)
db46 <- db46 %>%
  mutate(
    grupo3 = case_when(
      grupo2 %in% c("RE", "RE or SV", "SV") ~ "3. Dob_expo",
      grupo2 %in% c("SS1", "SS2", "V", "SS2 or V") ~ "2. Single_expo",
      grupo2 == "NVSN" ~ "1. NVSN",
      TRUE ~ NA_character_ # Para cualquier otro valor no especificado
    )
  )

m46 <- lm(volumen ~   cind_sexo+
            categoria_edad+grupo3 , data = db46)

tab_model(m46)

backwards = step(m46)

tab_model(backwards,show.aic = T)


############ Modelo L47
db47$razao_l45

m47a <- lm(volumen100 ~ cind_sexo+categoria_edad+
             razao_l45+
             razao_l46+ 
             factor(n.vacL47.x),
           data = db47)
summary(m47a)
tab_model(m47a,show.aic = T)



############ Modelo L48
db47$razao_l45
m48a <- lm(volumen100 ~ cind_sexo+categoria_edad+
             razao_l45+
             razao_l46+ 
             razao_l47+
             factor(n.vacL48.x)
           , data = db48)
tab_model(m48a,show.aic = T)


############ Modelo L49
db49
m49a <- lm(volumen100 ~ cind_sexo+
             categoria_edad+
             razao_l45+
             razao_l46+ 
             razao_l47+
             factor(n.vacL49.x)
           , data = db49)
tab_model(m49a,show.aic = T)





############ Modelo L45 

db45$slope


m45 <- lm(slope~ 
            cind_sexo +
            categoria_edad  , data = db45)

tab_model(m45,show.aic = T)


############ Modelo L46
db46$vac_cat.L46.x
m46 <- lm(slope ~  cind_sexo+
            categoria_edad+razao_l45+factor(vac_cat.L46.x), data = db46)

tab_model(m46)

library(dplyr)
db46 <- db46 %>%
  mutate(
    grupo3 = case_when(
      grupo2 %in% c("RE", "RE or SV", "SV") ~ "3. Dob_expo",
      grupo2 %in% c("SS1", "SS2", "V", "SS2 or V") ~ "2. Single_expo",
      grupo2 == "NVSN" ~ "1. NVSN",
      TRUE ~ NA_character_ # Para cualquier otro valor no especificado
    )
  )

m46 <- lm(volumen ~   cind_sexo+
            categoria_edad+grupo3 , data = db46)

tab_model(m46)

backwards = step(m46)

tab_model(backwards,show.aic = T)


############ Modelo L47
db47$razao_l45

m47a <- lm(slope ~ cind_sexo+categoria_edad+
             razao_l45+
             razao_l46+ 
             factor(n.vacL47.x),
           data = db47)
summary(m47a)
tab_model(m47a,show.aic = T)



############ Modelo L48
db47$razao_l45
m48a <- lm(slope ~ cind_sexo+categoria_edad+
             razao_l45+
             razao_l46+ 
             razao_l47+
             factor(n.vacL48.x)
           , data = db48)
tab_model(m48a,show.aic = T)


############ Modelo L49
db49
m49a <- lm(slope ~ cind_sexo+
             categoria_edad+
             razao_l45+
             razao_l46+ 
             razao_l47+
             factor(n.vacL49.x)
           , data = db49)
tab_model(m49a,show.aic = T)

