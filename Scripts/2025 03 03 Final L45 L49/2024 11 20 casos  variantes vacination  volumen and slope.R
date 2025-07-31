## Grafico inicial - Casos y variantes
## Juan P Aguilar Ticona 
## Version 2.0
setwd("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/99. Coorte de Pau da Lima/L50/Pau da Lima L50")

##########################################################################################
####### casos salvador 
##########################################################################################
##########################################################################################
####### F2A
##########################################################################################
library(readr)
library(dbplyr)
require(smooth)
require(Mcomp)

library(tidyverse)      # data manipulation and visualization
library(lubridate)      # easily work with dates and times
library(fpp2)           # working with time series data
library(zoo)            # working with time series data
library(scales)         # for function day_breaks
library(ggplot2)
library(dplyr)

a <- read_delim("Bancos/Covid curva Brasil/HIST_PAINEL_COVIDBR_2020_Parte1_28mai2024.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
b <- read_delim("Bancos/Covid curva Brasil/HIST_PAINEL_COVIDBR_2020_Parte2_28mai2024.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
c <- read_delim("Bancos/Covid curva Brasil/HIST_PAINEL_COVIDBR_2021_Parte1_28mai2024.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
d <- read_delim("Bancos/Covid curva Brasil/HIST_PAINEL_COVIDBR_2021_Parte2_28mai2024.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
e <- read_delim("Bancos/Covid curva Brasil/HIST_PAINEL_COVIDBR_2022_Parte1_28mai2024.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
f <- read_delim("Bancos/Covid curva Brasil/HIST_PAINEL_COVIDBR_2022_Parte2_28mai2024.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
h <- read_delim("Bancos/Covid curva Brasil/HIST_PAINEL_COVIDBR_2023_Parte1_28mai2024.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
i <- read_delim("Bancos/Covid curva Brasil/HIST_PAINEL_COVIDBR_2023_Parte2_28mai2024.csv", 
                delim = ";", escape_double = FALSE, trim_ws = TRUE)
j <- read_excel("Bancos/Covid curva Brasil/HIST_PAINEL_COVIDBR_2024_Parte1_14jun2024.xlsx")


# Banco do Brasil
Brasil <- rbind(a,b,c,d,e,f,h,i,j)

library(dplyr)
salvador <-  filter(Brasil, municipio == "Salvador")

salvador$DATA <- as.Date(salvador$data)
str(salvador$DATA)

max(salvador$DATA)
plot(salvador$DATA, salvador$casosNovos)
plot(salvador$data, salvador$obitosNovos)

salvador$CONFIRMADO <- salvador$casosNovos
salvador$OBITO <- salvador$obitosNovos


salvador$week <- strftime(salvador$DATA, format = "%U")
salvador$DATA <- as.POSIXct(salvador$DATA)
salvador$week <-format(salvador$DATA - as.difftime(weekdays(salvador$DATA), units = "days"), '%Y-%m-%d')


# Aggregate the date 

salvador$week <- floor_date(salvador$DATA, "week")

salvador2 <- salvador %>%                         # Aggregate data
  group_by(week) %>% 
  dplyr::summarize(cases = sum(casosNovos)) %>% 
  as.data.frame()

salvador3 <- salvador %>%                         # Aggregate data
  group_by(week) %>% 
  dplyr::summarize(deaths = sum(obitosNovos)) %>% 
  as.data.frame()


salvador_week <- merge(salvador2,salvador3,by.x = "week",all = T)
head(salvador_week)

salvador_week$week <- as.Date(salvador_week$week)
F2A <- ggplot() +
  geom_bar(data = salvador_week, aes(x=week,y=cases),stat="identity",fill = "#e9c46a", alpha = 0.75)+ theme_bw()+
  xlab("Date")+
  ylab("Number of cases")+
  scale_y_continuous(breaks = seq(0, 10000, by = 3000))+
  scale_x_date(limits = as.Date(c('2020-03-01','2024-04-01')),breaks=date_breaks("3 month"), labels=date_format("%b %y"))

F2A <- F2A + geom_bar(data = salvador_week, aes(x=week,y=deaths*8), stat = "identity", position = "stack",fill = "#dc0000ff", alpha=0.4) +
  scale_y_continuous(breaks = seq(0, 10000, by = 3000),sec.axis = dup_axis(trans = function(x) x/8, name = "Number of death", 
                                                                           breaks = seq(0,400,200)
  ))+
  theme_bw()  + theme_classic()+
  theme(axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),         
        axis.title.x=element_blank(),
        axis.title.y.left = element_text(color = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0, face = "plain"),
        axis.title.y.right = element_text(color = "grey20", size = 12, angle = 270, hjust = 0.5, vjust = 1.5, face = "plain"),
        #axis.text.x = element_blank(),
        axis.ticks.x= element_line(size = 1, linetype="dashed"),
        plot.margin=unit(c(0,0,0,0.5),"cm"))+ 
  geom_vline(xintercept = as.numeric(ymd("2019-09-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2019-11-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2020-11-18")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2021-02-26")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2021-07-14")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2021-10-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2022-03-25")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2022-08-31")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2022-11-16")), linetype="dashed", 
             color = "#82b0d2", size=1)+
 geom_vline(xintercept = as.numeric(ymd("2023-05-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2023-10-25")), linetype="longdash", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2024-03-24")), linetype="longdash", 
             color = "#82b0d2", size=1)
F2A


############################################
########### Variantes
############################################
### 
library(readr)
library(dplyr)
library(ggplot2)
library(lubridate)
library(scales)
library(tidyr)
library(readxl)

variantes <- read_excel("Bancos/1714050858125.metadata.tsv")

## Select salvador sequence
db <- subset(df, location == "Salvador" | location == "SALVADOR")

db1 <- db[,c("date","pangolin_lineage")]

table(db1$pangolin_lineage)
library(stringr)
str(db1$pangolin_lineage)


Ancenstral <- c("B.1", "B.1.1.28", "B.1.1.33", "B.1.356", "B.1.1", "B.1.1.378", "N.1", "P.2","N.9","P.7")
others <- c("Unassigned",
            "XDR","HV.1","HA.1","XBL.3","FL.4","FH.1","DS.2","CH.1.1","HK.3","JD.2")
db1$pangolin_lineage <- ifelse(db1$pangolin_lineage %in% Ancenstral, "Ancestral", db1$pangolin_lineage)
db1$pangolin_lineage <- ifelse(db1$pangolin_lineage %in% others, "Others", db1$pangolin_lineage)
db1$pangolin_lineage <- ifelse(db1$pangolin_lineage %in% c("EG.1","EG.6.1.1"), "EG.1 or EG.6.1.1", db1$pangolin_lineage)


db1$pangolin_lineage <- gsub("^BA\\.1.*", "BA.1*", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^BA\\.2.*", "BA.2*", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^BA\\.4.*", "BA.4*", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^BA\\.5.*", "BA.5*", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^BQ\\.1.1.*", "BQ.1.1*", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^XBB\\.1.5.*", "XBB.1.5*", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^XBB\\.1.18.*", "XBB.1.18*", db1$pangolin_lineage)

db1$pangolin_lineage <- gsub("^JD\\.1.*", "JD.1*", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^JN\\.1.*", "JN.1*", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^FE\\.1.*", "FE.1*", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^GK\\.1.*", "GK.1*", db1$pangolin_lineage)

db1$pangolin_lineage <- gsub("^AY\\..*", "Delta", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^BF\\..*", "Others", db1$pangolin_lineage)
db1$pangolin_lineage <- gsub("^BE\\..*", "Others", db1$pangolin_lineage)

db1$pangolin_lineage <- gsub("^P\\.1.*", "Gamma", db1$pangolin_lineage)

db1$pangolin_lineage <- ifelse(db1$pangolin_lineage == "B.1.617.2", "Delta", db1$pangolin_lineage)
#db1$pangolin_lineage <- ifelse(db1$pangolin_lineage == "P.1", "Gamma", db1$pangolin_lineage)
db1$pangolin_lineage <- ifelse(db1$pangolin_lineage == "B.1.1.7", "Alpha", db1$pangolin_lineage)
head(db1$pangolin_lineage)



## Wide dataframe
tabla_frecuencia <- table(db1$date, db1$pangolin_lineage)
db2 <- as.data.frame.matrix(tabla_frecuencia)
db2$date <- rownames(db2)
db2 <- db2[, c("date", colnames(db2)[-ncol(db2)])]

print(db2)

## month 
db2$date <- as.Date(db2$date)
db2$month <- as.Date(cut(db2$date, "month"))
db3 <- aggregate(. ~ month, db2, sum)
db3$date <- NULL
colnames(db3)[colnames(db3) == "JN.1*"] <- "JN.1"
db3$month <- as.Date(db3$month)
db3$JN.1


# Crear una nueva fila con valores NA y asignar valores específicos
new_row <- db3[1, ]  # Copia la estructura de db3
new_row[] <- NA       # Asigna NA a todas las columnas

# Asignar valores a las columnas específicas
new_row$month <- as.Date("2024-03-01")
new_row$JN.1 <- 30


# Agregar la nueva fila a db3
db3 <- bind_rows(db3, new_row)


# Crear una nueva fila con valores NA y asignar valores específicos
new_row <- db3[1, ]  # Copia la estructura de db3
new_row[] <- NA       # Asigna NA a todas las columnas

# Asignar valores a las columnas específicas
new_row$month <- as.Date("2024-04-01")
new_row$JN.1 <- 20
new_row$Others <- 1

# Agregar la nueva fila a db3
db3 <- bind_rows(db3, new_row)

db3[is.na(db3)] <- 0
####

library(tidyr)
variantes <- db3 %>% pivot_longer(!month,
                                  names_to='variant',
                                  values_to='Freq')

str(variantes$variant) 


library(dplyr)
variantes <- variantes %>% group_by(month) %>%
  mutate(percentage = Freq/sum(Freq)*100)

## orde de aparicion
variantes2 <- subset(variantes, Freq != 0)

# Obtener los niveles únicos de la variable 'variant' en el orden original
niveles_originales <- unique(variantes2$variant)
variantes2$variant <- factor(variantes2$variant, levels = niveles_originales)
order_var <- levels(variantes2$variant)
order_var
variantes$variant <- factor(variantes$variant, levels = order_var)
variantes$variant <- relevel(variantes$variant, "Others")



####  grafico 
library(viridis)
library(RColorBrewer)
library(ggsci)
library(viridis)
library(scico)
library(pals)

table(variantes$variant)


Fig2 <- ggplot(variantes, aes(x = month, y = percentage, fill = variant)) +
  geom_area()  +
  scale_x_date(limits = as.Date(c('2020-03-01','2024-03-30')),
               breaks=date_breaks("3 month"),
               labels=date_format("%b %y"))+
  theme_bw() +
  theme_classic() +
  labs(x = "", y = "Variants (%)") +
  theme(axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.title.x = element_blank(),
        axis.title.y.left = element_text(color = "grey20", size = 12, angle = 90, hjust = 0.5, vjust = 0, face = "plain"),
        axis.title.y.right = element_text(color = "grey20", size = 12, angle = 270, hjust = 0.5, vjust = 1.5, face = "plain"),
        axis.ticks.x = element_line(size = 1, linetype = "dashed"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "none",
        plot.margin = unit(c(0, 0, 0, 0.5), "cm"))+ 
  scale_fill_manual(values=as.vector(polychrome(27)))


library(ggplot2)
library(scales)

# Definir colores pastel para cada variante
colores_pastel2 <- c(
  "Others" = "#FFB3BA",
  "Ancestral" = "#FFDFBA",
  "Gamma" = "#FFFFBA",
  "Alpha" = "#BAFFC9",
  "Delta" = "#BAE1FF",
  "BA.1*" = "#D7B3FF",
  "BA.2*" = "#FFB3E6",
  "BA.4*" = "#B3E6FF",
  "BA.5*" = "#FFC3E1",
  "BQ.1.1*" = "#C3FF93",
  "BQ.1" = "#FF9E9D",
  "BQ.1.22" = "#FFC69D",
  "DL.1" = "#FFD89D",
  "XBB" = "#FFE39D",
  "XBB.1" = "#E3FF9D",
  "XBB.1.18*" = "#C5FF9D",
  "FE.1*" = "#9DFFB3",
  "XBB.1.15.1" = "#9DFFC5",
  "XBB.1.5*" = "#9DFFE3",
  "XBB.1.9" = "#9DFFFF",
  "EG.1 or EG.6.1.1" = "#9DE3FF",
  "GK.1*" = "#9DC5FF",
  "GK.3" = "#9DAFFF",
  "JD.1*" = "#B39DFF",
  "FL.1.5.1" = "#D19DFF",
  "XBB.1.16.6" = "#FF9DFF",
  "JN.1" = "#FF9DC5"
)

# Crear la figura con los colores asignados
Fig2 <- ggplot(variantes, aes(x = month, y = percentage, fill = variant)) +
  geom_area() +
  scale_x_date(limits = as.Date(c('2020-03-01', '2024-03-30')),
               breaks = date_breaks("3 month"),
               labels = date_format("%b %y")) +
  theme_bw() +
  theme_classic() +
  labs(x = "", y = "Variants (%)") +
  theme(axis.text.y = element_text(color = "grey20", size = 12),
        axis.text.x = element_text(angle = 30, hjust = 1, size = 12),
        axis.title.x = element_blank(),
        axis.title.y.left = element_text(color = "grey20", size = 12, angle = 90),
        axis.title.y.right = element_text(color = "grey20", size = 12, angle = 270),
        axis.ticks.x = element_line(size = 1, linetype = "dashed"),
        legend.title = element_blank(),
        legend.text = element_text(size = 10),
        legend.position = "top",
        plot.margin = unit(c(0, 0, 0, 0.5), "cm")) +
  scale_fill_manual(values = colores_pastel)

# Mostrar la figura
print(Fig2)




Fig2 <- Fig2+
  geom_vline(xintercept = as.numeric(ymd("2019-09-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2019-11-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2020-11-18")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2021-02-26")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2021-07-14")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2021-10-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2022-03-25")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2022-08-31")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2022-11-16")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2023-05-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2023-10-25")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2024-03-24")), linetype="dashed", 
             color = "#82b0d2", size=1)

Fig2

library(cowplot)

plot_grid(F2A, Fig2,  labels = c(), label_size = 12, ncol = 1, align = "vh", axis = "bt")



#plot_grid(F2A,
 #         Fig2, label_size = 12, ncol = 1, align = "vh", axis = "bt", labels = c("A","B","C","D","E"))
#PDF: 12*9 Landscape  



############################################
########### landscape
############################################

# Cargar las librerías necesarias
library(ggplot2)
library(scales)
library(lubridate)

# Datos
data <- data.frame(
  fecha = as.Date(c("2021-01-07", "2021-08-26", "2022-06-12", "2023-02-11")),
  volumen = c(7.7, 22.5, 58.9, 64.3),
  slope = c(0.31, 0.41, 0.47, 0.44)
)

# Fechas de las líneas verticales
vertical_lines <- as.Date(c("2019-09-09", "2019-11-09", "2020-11-18", "2021-02-26", 
                            "2021-07-14", "2021-10-09", "2022-03-25", "2022-08-31", 
                            "2022-11-16", "2023-05-09", "2023-10-25","2024-03-24"))

# Crear la figura con colores pastel y líneas verticales
Fig3 <- ggplot(data) +
  # Eje para el volumen con color pastel
  geom_line(aes(x = fecha, y = volumen), color = "#A3C4BC", size = 1) +  # Pastel verde-azulado
  geom_point(aes(x = fecha, y = volumen), color = "#A3C4BC", size = 3) +
  # Eje para el slope con color pastel
  geom_line(aes(x = fecha, y = slope * 100), color = "#F7B7A3", size = 1) +  # Pastel rosado
  geom_point(aes(x = fecha, y = slope * 100), color = "#F7B7A3", size = 3) + 
  scale_x_date(date_breaks = "3 months", date_labels = "%m-%Y",
               limits = as.Date(c('2020-03-01', '2023-06-01')))+
  # Añadir líneas verticales
  lapply(vertical_lines, function(x) {
    geom_vline(xintercept = as.numeric(x), linetype = "dashed", color = "#82b0d2", size = 1)
  }) +
  scale_y_continuous(
    name = "Volumen (%)",
    sec.axis = sec_axis(~ . / 100, name = "Slope")
  ) +
  # Etiquetas y título
  labs(
    x = "",
    title = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )+theme_bw()+theme_classic()


plot_grid(F2A,
          Fig2, Fig3, label_size = 12, ncol = 1, align = "vh", axis = "bt", labels = c("A","B","C","D","E"))



############
##### vacination 
#############

library(readxl)
vac <- read_excel("C:/Users/Pablo/OneDrive - FIOCRUZ/Pablo - Fiocruz/1. MSD/1. Landscape clean/MSD/dataset/bases anteriores/2024 06 05 Vacina covid resumen.xlsx")


library(ggplot2)
library(dplyr)

# Asegurarte de que las fechas están en formato Date
vac$data1 <- as.Date(vac$data1)
vac$data2 <- as.Date(vac$data2)
vac$data3 <- as.Date(vac$data3)
vac$data4 <- as.Date(vac$data4)

# Calcular el porcentaje acumulado para el primer conjunto de datos
data_cum1 <- data.frame(date = vac$data1) %>%
  count(date) %>% 
  arrange(date) %>%
  mutate(
    cumulative = cumsum(n),
    percentage = (cumulative / 613) * 100 # Ajustar el total
  )

# Calcular el porcentaje acumulado para el segundo conjunto de datos
data_cum2 <- data.frame(date = vac$data2) %>%
  count(date) %>% 
  arrange(date) %>%
  mutate(
    cumulative = cumsum(n),
    percentage = (cumulative / 613) * 100 # Ajustar el total
  )
# Calcular el porcentaje acumulado para el segundo conjunto de datos
data_cum3 <- data.frame(date = vac$data3) %>%
  count(date) %>% 
  arrange(date) %>%
  mutate(
    cumulative = cumsum(n),
    percentage = (cumulative / 613) * 100 # Ajustar el total
  )
# Calcular el porcentaje acumulado para el segundo conjunto de datos
data_cum4 <- data.frame(date = vac$data4) %>%
  count(date) %>% 
  arrange(date) %>%
  mutate(
    cumulative = cumsum(n),
    percentage = (cumulative / 613) * 100 # Ajustar el total
  )






vac_combined <- ggplot() +
  # Primera curva
  # Primera curva
  geom_area(data = data_cum1, aes(x = date, y = percentage), fill = "skyblue", alpha = 1) +
  geom_line(data = data_cum1, aes(x = date, y = percentage), color = "skyblue", size = 1.2) +
  # Segunda curva
  geom_area(data = data_cum2, aes(x = date, y = percentage), fill = "pink", alpha = 1) +
  geom_line(data = data_cum2, aes(x = date, y = percentage), color = "pink", size = 1.2) +
  # Segunda curva
  geom_area(data = data_cum3, aes(x = date, y = percentage), fill = "#107F97", alpha = 1) +
  geom_line(data = data_cum3, aes(x = date, y = percentage), color = "#107F97", size = 1.2) +
  # Segunda curva
  geom_area(data = data_cum4, aes(x = date, y = percentage), fill = "#BDCE45", alpha = 1) +
  geom_line(data = data_cum4, aes(x = date, y = percentage), color = "#BDCE45", size = 1.2) +
  #Etiquetas y escalas
  labs(
    title = "Curvas Acumulativas de Vacunación",
    x = "Fecha de Vacunación",
    y = "Porcentaje Acumulado"
  ) +
  scale_x_date(
    date_breaks = "3 months", 
    date_labels = "%m-%Y"
  ) +
  scale_y_continuous(
    limits = c(0, 100), 
    labels = scales::percent_format(scale = 1)
  ) +
  scale_x_date(limits = as.Date(c('2020-03-01', '2024-03-30')),
               breaks = date_breaks("3 month"),
               labels = date_format("%b %y")) + # Corta visualmente los datos
  theme_classic()+ 
  geom_vline(xintercept = as.numeric(ymd("2019-09-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2019-11-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2020-11-18")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2021-02-26")), linetype="dashed", 
             color = "#82b0d2", size=1)+ 
  geom_vline(xintercept = as.numeric(ymd("2021-07-14")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2021-10-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2022-03-25")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2022-08-31")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2022-11-16")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2023-05-09")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2023-10-25")), linetype="dashed", 
             color = "#82b0d2", size=1)+
  geom_vline(xintercept = as.numeric(ymd("2024-03-24")), linetype="dashed", 
             color = "#82b0d2", size=1)

vac_combined


plot_grid(F2A,
          Fig2, vac_combined, label_size = 12, ncol = 1, align = "vh", axis = "bt", labels = c("A","B","C","D","E"))

