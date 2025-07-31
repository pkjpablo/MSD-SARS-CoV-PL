library(writexl)
write_xlsx(db2, "dataset/bancos general/192 Pau da Lima all.xlsx")


## primera parte 
## 

#Clear existing data and graphics
rm(list=ls())
graphics.off()
#Longitudinal samples

library(readr)
library(readxl)
library(ggsci)
library(plyr)
library(table1)
library(writexl)


# abrir y juntar bancos 
db <- read_csv("dataset/MSDresults_Salvador_031824.csv")
db32new <- read_excel("dataset/MSDresults_Salvador_051524 new results cate.xlsx", 
                      sheet = "Panel 32")
db23new <- read_excel("dataset/MSDresults_Salvador_051524 new results cate.xlsx", 
                      sheet = "Panel 23")
db23new$panel[db23new$panel == "Panel 24"] <- "Panel 23"
db23new <- subset(db23new,!plate=="Plate 23")


table(db32new$survey)
table(db23new$survey)

a <- rbind.fill(db, db23new)
db <- rbind.fill(a, db32new)

table(db$dilution)
table(db$panel)

# sera ace2 Pl

db <- subset(db, survey != "VIGIPL")
db <- subset(db, sample_type == "Serum")
db <- subset(db, test_type == "ACE2")
db <- subset(db, dilution == "1:20")
db <- subset(db, survey != "L38")
db <- subset(db, is.na(cohort) | cohort != "comparar 23")

table(db$cohort)
table(db$plate)
table(db$panel)
table(db$survey,db$panel)


### create a unica id 
db$idjp <- paste(db$survey, db$indid,  db$sample_type, sep = "_")

table(db$survey)
table(db$panel)
table(db$survey,db$panel)

###  filtrar panel 23
p23 <- subset(db, panel == "Panel 23" | panel == "Panel 24")
table(p23$survey,p23$panel)

p23$panel[p23$panel == "Panel 24"] <- "Panel 23"
table(p23$survey,p23$panel)

# Seleccionar todas las columnas desde p23$pct.ancestral_p23 hasta p23$pct.P.1 y la columna p23$idjp
p23 <- p23[, c((which(names(p23) == "idjp")),
               (which(names(p23) == "pct.ancestral_p23")):which(names(p23) == "pct.P.1"))]
head(p23)


# Identificar las filas duplicadas en la variable p23$odjp
duplicados_idjp <- duplicated(p23$idjp) | duplicated(p23$idjp, fromLast = TRUE)
# Crear un subset con las filas duplicadas
subset_duplicados <- subset(p23, duplicados_idjp)

library(dplyr)

p23mean <- aggregate(. ~ idjp, data = p23, FUN = function(x) mean(x, na.rm = TRUE))

p23_unique <- p23[!duplicated(p23$idjp, fromLast = TRUE), ]
p23mean <- p23_unique

p23mean$survey <- substr(p23mean$idjp, 1, 3)
p23mean$idjp <- substr(p23mean$idjp, 5, nchar(p23mean$idjp))




library(tidyr)
p23 <- p23mean %>%
  pivot_wider(
    names_from = survey,
    names_sep = ".",
    values_from = c(pct.ancestral_p23:pct.P.1)
  )

#p23L45 <- subset(p23mean,survey=="L45")
#p23L46 <- subset(p23mean,survey=="L46")

table(p23mean$survey)
###################### panel 32

#Longitudinal samples

library(readr)
db <- read_csv("dataset/MSDresults_Salvador_031824.csv")
db32new <- read_excel("dataset/MSDresults_Salvador_051524 new results cate.xlsx", 
                      sheet = "Panel 32")
db23new <- read_excel("dataset/MSDresults_Salvador_051524 new results cate.xlsx", 
                      sheet = "Panel 23")

a <- rbind.fill(db, db23new)
db <- rbind.fill(a, db32new)

## Cuantas amostra de soro de L45 46 p23 e L47 L48 p32 forma testadas 

db <- subset(db, survey != "VIGIPL")
db <- subset(db, sample_type == "Serum")
db <- subset(db, test_type == "ACE2")
db <- subset(db, dilution == "1:20")
db <- subset(db, survey != "L38")


### create a unica id 

db$idjp <- paste(db$survey, db$indid,  db$sample_type, sep = "_")
p32 <- subset(db, panel == "Panel 32")
table(p32$survey)

# Seleccionar todas las columnas desde p32$pct.ancestral_p32 hasta p32$pct.P.1 y la columna p32$idjp
p32 <- p32[, c( (which(names(p32) == "idjp")),
                (which(names(p32) == "pct.ancestral_p32")):(which(names(p32) == "pct.XBB.1")))]

head(p32)

# Identificar las filas duplicadas en la variable p32$odjp
duplicados_idjp <- duplicated(p32$idjp) | duplicated(p32$idjp, fromLast = TRUE)

# Crear un subset con las filas duplicadas
subset_duplicados <- subset(p32, duplicados_idjp)

library(dplyr)

p32mean <- aggregate(. ~ idjp, data = p32, FUN = mean, na.rm = TRUE)
p32mean$survey <- substr(p32mean$idjp, 1, 3)
p32mean$idjp <- substr(p32mean$idjp, 5, nchar(p32mean$idjp))


library(tidyr)
p32 <- p32mean %>%
    pivot_wider(
    names_from = survey,
    names_sep = ".",
    values_from = c(pct.ancestral_p32:pct.XBB.1)
  )

all <- merge(p23,p32,by = "idjp",all = T)

# Cargar el paquete dplyr
library(dplyr)

# Contar las filas con algún valor en las columnas que contienen "ancestral" en su nombre
all <- all %>%
  mutate(number = rowSums(!is.na(select(., contains("ancestral")))))

all <- all %>%
  mutate(number23 = rowSums(!is.na(select(., contains("ancestral_p23")))))
all <- all %>%
  mutate(number32 = rowSums(!is.na(select(., contains("ancestral_p32")))))

table1(~factor(number23)+factor(number32)|factor(number),all)

all <- all %>%
  mutate(idnova = substr(idjp, 1, 9))

library(writexl)

write_xlsx(all, 'dataset/2025 02 21 All ACE2 resultsv2.xlsx')

#write_xlsx(all, 'dataset/2024 06 17 All ACE2 resultsv2.xlsx')

#write_xlsx(all, 'dataset/2024 06 17 All ACE2 results.xlsx')

#### Salvando bancos 

# Para landscape study 
table(all$number)
dbland <- subset(all,number=="10")

table(all$number)
vars_to_exclude <- c("number", "number23","number32")
dbland <- dbland[, setdiff(names(dbland), vars_to_exclude)]
last_col <- names(dbland)[ncol(dbland)]
new_order <- c(names(dbland)[1], last_col, names(dbland)[-c(1, ncol(dbland))])
dbland <- dbland[, new_order]


#write_xlsx(dbland, 'dataset/2024 06 05 137 participantsv2.xlsx')

#write_xlsx(dbland, 'dataset/2025 01 07 205 participantsv2.xlsx')
write_xlsx(dbland, 'dataset/2025 02 25 205 participantsv2.xlsx')


## Sociodemografico
library(readxl)
L45socio <- read_excel("dataset/bases anteriores/L45 sociodemografico pablo.xlsx")

socio <- L45socio %>%
  select(idnova,cens_idade, idade_cat, cind_sexo, raca2, serie, schooling, estado_civil3, trabalho2)

all2 <- merge(dbland,socio,by = "idnova",all.x = T)

#write_xlsx(all2, 'dataset/2024 06 05 204 participantsv2.xlsx')

library(gtsummary)
all2
tab_socio <- all2 %>% 
  select(cens_idade,idade_cat,cind_sexo) %>% 
  tbl_summary()

tab_socio

## 
library(readxl)

iggL45 <- read_csv("dataset/bases anteriores/L45SoroinqueritoHarm-L45Basico_DATA_2024-03-07_1402.csv")
iggL45 <- subset(iggL45,!is.na(razao_l45))
iggL45 <- subset(iggL45,!is.na(idnova))
iggL45 <- iggL45 %>%
  select(idnova,dtcoleta,razao_l45)

iggL46 <- read_csv("dataset/bases anteriores/L46IFAASHarmonizado-L46Basico_DATA_2024-03-07_1403.csv")
iggL46 <- subset(iggL46,!is.na(razao_l46))
iggL46 <- subset(iggL46,!is.na(idnova))
iggL46 <- iggL46 %>%
  select(idnova,dtcoleta,razao_l46)

igg <- merge(iggL45,iggL46,by = "idnova",suffixes = c("L45","L46"), all = T)

all2 <- merge(all2,igg,by = "idnova",all.x = T)

#all3 <- subset(all2,number=="8")


tab_socio <- all2%>% 
  select(cens_idade,idade_cat,cind_sexo) %>% 
  tbl_summary()

tab_socio


## PRNT 
prnt <- read_excel("dataset/SalvadorL45L46Neut022423 (3).xlsx")

prnt <- prnt[,c("indid","PRNT50.Ancestral.L45",
                "PRNT50.Gamma.L45",
                "PRNT50.Delta.L45",
                "PRNT50.Omicron.L45",
                "PRNT50.Ancestral.L46",
                "PRNT50.Gamma.L46",
                "PRNT50.Delta.L46",
                "PRNT50.Omicron.L46")]


prnt2 <- prnt %>%
  filter(complete.cases(.))

all3 <- merge(all2,prnt2, by.x = "idnova",by.y = "indid",all.x = T)


## Vaccine
library(readxl)
vac <- read_excel("dataset/bases anteriores/draft/2024 05 10 Vacina covid.xlsx", 
                                       col_types = c("numeric", "text", "text", 
                                                     "text", "numeric", "date", "text", 
                                                     "text", "text", "numeric", "date", 
                                                     "text", "date", "text", "date", "text", 
                                                     "date", "text", "date", "text", "date", 
                                                     "text"))

vac <- vac[, c(1, 2, 6, 7, 9:22)]

all4 <- merge(all3,vac, by.x = "idnova",by.y = "idnova",all.x = T)


library(readxl)
expo <- read_excel("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/Antigenic maps/R antigenic maps/Bancos/exposure.xlsx")
expo$idnova
all4<- merge(all4,expo,by="idnova",all.x = T)


###################
# Estudo de Emilia
###################

p23L45 <- subset(p23mean,survey=="L45")
p23L46 <- subset(p23mean,survey=="L46")

p32L45 <- subset(p32mean,survey=="L45")
p32L46 <- subset(p32mean,survey=="L46")


emilia <- subset(all4, !is.na(PRNT50.Ancestral.L45))
emilia <- subset(emilia, !is.na(pct.ancestral_p23.L45))
emilia <- subset(emilia, !is.na(pct.ancestral_p32.L46))
emilia <- subset(emilia, !is.na(pct.ancestral_p32.L46))

table(emilia$conf,useNA = "always")

library(table1)
limites_edad <- c(-Inf,18, 29, 44, 59, Inf)
etiquetas_edad <- c("<18","18-29", "30-44", "45-59", "60+")
emilia$categoria_edad <- cut(emilia$cens_idade, breaks = limites_edad, labels = etiquetas_edad, right = FALSE)

emilia$razao_l45

# Definir las variables
variables <- c("razao_l45", "razao_l46")
# Iterar sobre cada variable
for (var in variables) {
  # Crear nueva variable categórica
  emilia[[paste0(var, "_grupo")]] <- ifelse(emilia[[var]] > 0.8, "Positive", "Negative")
}


# Reemplazar los valores en emilia$vac_covid
emilia$razao_l45_grupo <- ifelse(emilia$razao_l45_grupo == "Positive",
                             "SS1",
                             ifelse(emilia$razao_l45_grupo == "Negative",
                                    "NoSS1", emilia$razao_l45_grupo))


# Reemplazar los valores en emilia$vac_covid
emilia$razao_l46_grupo <- ifelse(emilia$razao_l46_grupo == "Positive",
                             "SS2",
                             ifelse(emilia$razao_l46_grupo == "Negative",
                                    "NoSS2", emilia$razao_l46_grupo))

#### vac
# Calcular la diferencia en días entre las fechas
emilia$vacL46 <- difftime(emilia$data1, emilia$dtcoletaL46, units = "days")

# Marcar como "not vaccinated" si alguna fecha es NA
emilia$vacL46[is.na(emilia$data1)] <- "not vaccinated"
emilia$vacL46[emilia$vacL46 > 0] <- "not vaccinated"
emilia$vacL46[emilia$vacL46 < 0.1] <- "vaccinated"

#### reinfectado      ----

emilia$reinfectado <- emilia$razao_l46/emilia$razao_l45
hist(emilia$reinfectado)

# Definir el punto de corte
punto_de_corte <- 1.5

# Crear la nueva variable dn$RE
emilia$RE <- ifelse(emilia$reinfectado >= 1.5 & emilia$razao_l45_grupo == "SS1", "RE", "NoRE")
table(emilia$RE)


emilia$grupo <- paste(emilia$razao_l45_grupo,
                      emilia$razao_l46_grupo,
                      emilia$RE,
                      emilia$vacL46)

table1(~grupo, emilia)

table1(~categoria_edad+cind_sexo+raca2+schooling+
         trabalho2+razao_l45_grupo+razao_l46_grupo+
         razao_l46_grupo+RE+vacL46+exposure_group|exposure_group,emilia)


table1(~categoria_edad+cind_sexo+raca2+schooling+
         trabalho2+razao_l45_grupo+razao_l46_grupo+
         razao_l46_grupo+exposure_group+RE+vacL46,emilia)




a <- emilia$idnova

saveRDS(a,file = "dataset/51 emilia")

