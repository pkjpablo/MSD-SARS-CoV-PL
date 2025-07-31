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




df45$volumen <- df45$volumen/270*100
df46$volumen <- df46$volumen/270*100
df47$volumen <- df47$volumen/270*100
df48$volumen <- df48$volumen/270*100
df49$volumen <- df49$volumen/270*100

df45$id <- substr(df45$id, 5, nchar(df45$id))
df46$id <- substr(df46$id, 5, nchar(df46$id))
df47$id <- substr(df47$id, 5, nchar(df47$id))
df48$id <- substr(df48$id, 5, nchar(df48$id))
df49$id <- substr(df49$id, 5, nchar(df49$id))


a <- merge(df45, df46, by = "id", suffixes = c("L45", "L46"))
a$dif <- a$volumenL46-a$volumenL45
a$period <- "L45-L46"

b <- merge(df46, df47, by = "id", suffixes = c("L46", "L47"))
b$dif <- b$volumenL47-b$volumenL46
b$period <- "L47-L46"

d <- merge(df47, df48, by = "id", suffixes = c("L47", "L48"))
d$dif <- d$volumenL48-d$volumenL47
d$period <- "L48-L47"

e <- merge(df48, df49, by = "id", suffixes = c("L48", "L49"))
e$dif <- e$volumenL49-e$volumenL48
e$period <- "L49-L48"


library(plyr)
dbdifv <- rbind.fill(a, b, d, e)

dbx <- read_excel("dataset/bases anteriores/socio2.xlsx")



library(ggplot2)

ggplot(dbdifv, aes(x = factor(period), y = dif, color = factor(period), fill = factor(period))) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Oculta los outliers para mejorar visualización
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +  # Puntos dispersos
  labs(x = "", y = "Δ % Inhibition", title = "") +
  theme_classic() +  # Tema clásico
  theme(axis.text.x = element_text(angle = 60, hjust = 1))

library(ggplot2)

ggplot(dbdifv, aes(x = factor(period), y = dif, color = factor(period), fill = factor(period))) +
  geom_boxplot(alpha = 0.3, outlier.shape = NA) +  # Oculta los outliers
  geom_jitter(width = 0.2, size = 1, alpha = 0.5) +  # Puntos dispersos
  labs(x = "Period", y = "Δ % Inhibition", title = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  scale_color_manual(values = c("#2b6caa", "#a98122", "#1b5467", "#4c4471")) +  # Bordes
  scale_fill_manual(values = c("#2b6caa", "#a98122", "#1b5467", "#4c4471"))  # Rellenos
