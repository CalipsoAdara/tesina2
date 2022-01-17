# MJO con regiones
# Analizar si las inicializaciones con MJO activos afecta a los scores segun la$
# Analisis zona SACZ y SUR PATAGONIA
# # Si resto activo - inactivo
# donde rmse sea positivo siginifica que fue mayor en activo ----> MAL
# donde acc sea positivo significa que fue mayor en activo -----> BIEN
# 
#
# by Lucia M Castro
#------------------------------------------------------------------------------$

# Limpio enviroment
rm(list=ls())

# Seteo directorio
svpath = "/home/lucia.castro/tesina2/resultados"
setwd("/vegeta/datos/SubX/mjo.lucia.cas")

# Cargo paquetes
library(secr)
library(reshape2)
library(csv)
library(dplyr)
library(data.table)
library(ggplot2)
library(metR)



# Cargo mis funciones
source("/vegeta/datos/SubX/mjo.lucia.cas/funciones.R")

groups=c('GMAO','RSMAS','ESRL','ECCC','NRL','EMC','MME')                       
models=c('GEOS_V2p1','CCSM4','FIMr1p1','GEM','NESM','GEFS','SAT')   

# Cargo los datos de eventos
df_rmm <- readRDS("/vegeta/datos/SubX/mjo.lucia.cas/df_rmm.rds")
df_eventos <- readRDS("/vegeta/datos/SubX/mjo.lucia.cas/df_eventos.rds")
fechas_act <- as.character(df_rmm$DATE)

# Poligonos. Lon de menor a mayor, el primer punto se repite para cerrar el poligono
SP <- data.frame(x_coords = c(291,288,291,298,291),
                 y_coords = c(-30,-40,-53,-40,-30))

SACZ <- data.frame(x_coords = c(305,305,310,321,305),
                   y_coords = c(-10,-25,-30,-10,-10))

# Separar segun la fase inicial
# Bins = [8,1] [2,3] [4,5] [6,7]

# Debo calcular las metricas (RMSE, ACC) con solo los datos que caen en el bin de fase inicial. A estos luego se restan
# la metrica calculada para las fechas inactivas de MJO. En mjo_fases.R ya se calculo y guardaron los datos.

Bins = levels(df_eventos$Bin)
listagraficosRMSE <- list()
listagraficosACC <- list()


for (g in 1:length(groups)) {
  grupo = groups[g]
  model = models[g]
  
  metricINA <- readRDS(paste0("./metricsMJO_inact_",grupo,".rds"))[c(1,3)] #Leo rmse y acc solo
  
  for (b in Bins) { # por cada Bin
    # Cargo los datos
    metricACT <- readRDS(paste0("./ScoresBins/",grupo,b))
    # Resto
    resta_bin <- Map('-', metricACT, metricINA) #Resto ambas listas 
    
    for(m in length(resta_bin)) { # Por cada metrica (acc y rmse)
      # Convierto a data frame
      df_resta <- reshape2::melt(resta_bin[[m]])
      
      # Restringir el data frame al area del poligono (primeras 2 col son lat y lon)
      puntossacz=pointsInPolygon(df_resta[,1:2],SACZ) 
    }
  
    
    
  }# End Bin
}# End model