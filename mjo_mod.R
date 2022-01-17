# MJO analisis con los modelos. 

# ---------------------------------------------------------------------------------
# Limpiar el espacio
rm(list=ls())

# cargar paquetes
library(metR)
library("grid")
library(RNetCDF)
library(raster)
library('pracma')
library("gridExtra")
library(data.table)
library("ggplot2")


# Cargo mis funciones
source("/home/lucia.castro/tesina/funciones.R")

# Seteo el directorio
setwd("/home/lucia.castro/SubX_processed_Rdata")

groups=c('GMAO','RSMAS','ESRL','ECCC','NRL','EMC','MME')                       
models=c('GEOS_V2p1','CCSM4','FIMr1p1','GEM','NESM','GEFS','SAT')   

var_mjo = c("RMM1","RMM2","amplitude","phase")


# # DATOS DE RMM DE LOS MODELOS 
# # Directorio donde guardar
# outPath = '/datos/SubX/hindcast/rmm/'
# 
# inFname = "/datos/SubX/hindcast/rmm/Amplitude/NRL_amp.nc"
# data = metR::ReadNetCDF(inFname, out='array')
# data = data[[1]]
# 
# dimnames(data)$S <- as.character(seq.Date(as.Date("1999-01-01"),as.Date("2015-12-30"),1))


# Supongo que donde hubo MJO activo para las observaciones, el modelo tambien tiene evento activo
# Me interesa solo donde se inicio el modelo con evento activo
# Podria separar por fases al iniciar 

# Cargo los datos de eventos
df_rmm <- readRDS("/vegeta/datos/SubX/mjo.lucia.cas/df_rmm.rds")
df_eventos <- readRDS("/vegeta/datos/SubX/mjo.lucia.cas/df_eventos.rds")
fechas_act <- as.character(df_rmm$DATE)
# Cargo datos del rho1 para la significancia
rho1 <- readRDS("./rho1.rds")

for (g in 1:length(groups)) {
  grupo = groups[g]
  model = models[g]
  # Busco que startdates coinciden con los eventos activos
  TargetDate <- readRDS(paste0("./targetdate_",grupo,"_ONDEFM.rds"))
  startdate = dimnames(TargetDate)$startdate
  posMJO = startdate %in% fechas_act
  
  modelweek<-readRDS(paste0("./modelweek_",grupo,".rds"))
  obsweek<-readRDS(paste0("./obsweek_",grupo,".rds"))
  
  modMJO = modelweek[,,!posMJO,]
  obsMJO = obsweek[,,!posMJO,]
  
  # Cantidad de inicios antes y despues de restringir en los eventos
  nstartdate = length(startdate)
  nstartdateMJO = sum(posMJO)
  
  # Calculo las distintas métricas por cada lon/lat/targetweek
  dif = (obsMJO - modMJO)
  me = apply(dif, c(1,2,4),FUN = mean, na.rm = TRUE)
  mae = apply(abs(dif), c(1,2,4), FUN = mean, na.rm = TRUE) 
  rmse = sqrt(apply(dif^2,c(1,2,4), FUN = mean, na.rm = TRUE))
  desvio = apply(dif,c(1,2,4),FUN = sd, na.rm = TRUE)
  var = (1-sqrt(rmse))/desvio
  acc = ACC(Lon=66, Lat = 76, Model = modMJO, Anom = obsMJO) # tarda un cacho
  
  # Calculo estadistico de prueba
  # RHO1 sirve para todas las weeks (1,2.3 y 4). Repito el tamaño de muestra para cada semana
  n_eff = nstartdate*((1 - rho1)/(1 + rho1))
  n_eff = array(n_eff, dim = c(66,76,4))
  t = (acc * sqrt(n_eff - 2)) / sqrt(1-acc^2)
  
  # Significancia de 0.05 y grados de libertad 
  critc = qt(p=0.95, df = trunc(n_eff))
  test = t < critc
  
  # Renombro dimensiones 
  dt.me = ggScoreSemanal(me)
  dt.mae = ggScoreSemanal(mae)
  dt.rmse = ggScoreSemanal(rmse)
  dt.var = ggScoreSemanal(var)
  dt.acc = ggScoreSemanal(acc)
  
  metrics <- list(rmse,me,acc,var)
  saveRDS(metrics, paste0("./metricsMJO_inact_",grupo,".rds"))
  #---------------------------------------------------------------------------------------
  #  Gráficos  
  #---------------------------------------------------------------------------------------
  g1 <- GraphDiscreteMultiple(Data = dt.rmse, Breaks = seq(0,3,0.25),Label = "RMSE",Paleta = "YlOrRd", Direccion = "1")
  g2 <- GraphDiscreteMultiple(Data = dt.me, Breaks = seq(-0.1,0.1,0.025), Label = "ME",Paleta = "RdBu",Direccion = "-1")
  g3 <- GraphMultiplePuntos(Data = dt.acc, ArLogic = test, Breaks = seq(0,1,0.20), Label = "ACC",Paleta = "YlGn",Direccion = "1")
  g4 <- GraphDiscreteMultiple(Data = dt.var, Breaks = seq(-0.5,0.5,0.10), Label = "NRMSE",Paleta = "RdBu",Direccion = "-1")
  
  
  
  title = title <- paste("SubX ",grupo,"-",model," Inicios",nstartdateMJO,"/",nstartdate,
                         "\nMJO inactive events tasa (99-15, Oct-Mar) ")
  fig <- grid.arrange(g1,g2,g3,g4, ncol = 1,top = textGrob(title,gp=gpar(fontsize=13,font=3)))
  ggsave(filename=paste0("/home/lucia.castro/SubX_processed_Rdata/scores_map_MJO_inac_",grupo,".png"),plot=fig,width = 10, height = 15)
  
}


# Hacer la resta entre los scores inactivos y cuando hay evento activo
# TOTAL

for (g in 1:length(groups)) {
  grupo = groups[g]
  model = models[g]
  # Busco que startdates coinciden con los eventos activos
  TargetDate <- readRDS(paste0("./targetdate_",grupo,"_ONDEFM.rds"))
  startdate = dimnames(TargetDate)$startdate
  posMJO = startdate %in% fechas_act
  # Cantidad de inicios antes y despues de restringir en los eventos
  nstartdate = length(startdate)
  nstartdateMJO = sum(posMJO)
  
  # Cargo los scores
  metricINA <- readRDS(paste0("./metricsMJO_inact_",grupo,".rds"))
  metricACT <- readRDS(paste0("./metricsMJO_",grupo,".rds"))
  
  # Hago la resta de cada metrica
  metricANOM <- list()
  for (m in 1:length(metricMOD)) {
    metric <- metricACT[[m]] - metricINA[[m]]
    df <- ggScoreSemanal(metric)
    metricANOM[[m]] <- df
  }
  #---------------------------------------------------------------------------------------
  #  Gráficos  
  #---------------------------------------------------------------------------------------
  
  test = array(NA, c(66,76,4))
  # Como es (todo - mjo) significa que cuando es negativo el mjo tiene valor mayor
  
  g1 <- GraphDiscreteMultiple(Data = metricANOM[[1]], Breaks = seq(-0.2,0.2,0.05),Label = "RMSE",Paleta = "RdBu", Direccion = "-1")
  g3 <- GraphMultiplePuntos(Data = metricANOM[[3]], ArLogic = test, Breaks = seq(-0.2,0.2,0.05), Label = "ACC",Paleta = "RdBu",Direccion = "-1")
  
  
  
  title = paste0("SubX ",grupo,"-",model,"\nMJO activo - inactivo (99-15, Oct-Mar)")
  fig <- grid.arrange(g1,g3, ncol = 1,top = textGrob(title,gp=gpar(fontsize=13,font=3)))
  ggsave(filename=paste0("/home/lucia.castro/SubX_processed_Rdata/MJO/scores_map_MJO_dff",grupo,".png"),plot=fig,width = 10, height = 15)
  
}

# Si resto activo - inactivo
# donde rmse sea positivo siginifica que fue mayor en activo ----> MAL
# donde acc sea positivo significa que fue mayor en activo -----> BIEN
# Seteo los parametros de mapa y gradiente 

# Separar segun la fase inicial
# Bins = [8,1] [2,3] [4,5] [6,7]
# Quiero hacer la diferencia entre activos - inactivos para C/bin

Bins = levels(df_eventos$Bin)
nstartdateMJOBIN <- list()
listagraficosRMSE <- list()
listagraficosACC <- list()


for (g in 1:length(groups)) {
  grupo = groups[g]
  model = models[g]
  
  # Cargo las fechas del modelo
  TargetDate <- readRDS(paste0("./targetdate_",grupo,"_ONDEFM.rds"))
  startdate = dimnames(TargetDate)$startdate
  # Cargo datos de modelo y observaciones
  modelweek<-readRDS(paste0("./modelweek_",grupo,".rds"))
  obsweek<-readRDS(paste0("./obsweek_",grupo,".rds"))
  metricINA <- readRDS(paste0("./metricsMJO_inact_",grupo,".rds"))[c(1,3)] #Leo rmse y acc solo
  
  for (b in Bins) { # por cada Bin
    # Busco que startdates coinciden con los eventos activos
    fechas_act_bin <- as.character(df_rmm[Bin==b,DATE])
    posMJOBIN = startdate %in% fechas_act_bin
    
    # Cantidad de inicios antes y despues de restringir en los eventos
    # informacion para luego crear tabla
    nstartdate = length(startdate)
    nstartdateMJOBIN <- c(nstartdateMJOBIN,sum(posMJOBIN))
    
    # Calculo las distintas métricas por cada lon/lat/targetweek
    sco_bin<-Metrics(obsweek[,,posMJOBIN,],modelweek[,,posMJOBIN,])
    saveRDS(sco_bin,paste0("./MJO/ScoresBins/",grupo,b))
    resta_bin <- Map('-', sco_bin, metricINA) #Resto ambas listas 
    
    # Graficos
    g1 <- GraphDiscreteMultiple(Data = ggScoreSemanal(resta_bin[[1]]), Breaks = seq(-0.2,0.2,0.05),Label = "RMSE",Paleta = "RdBu", Direccion = "1") 
    g2 <- GraphDiscreteMultiple(Data = ggScoreSemanal(resta_bin[[2]]), Breaks = seq(-0.2,0.2,0.05), Label = "ACC",Paleta = "RdBu",Direccion = "-1")
    
    # Lo completo asi en vez de append porque sino guarda listas en vez de ggplots
    len <- length(listagraficosRMSE)
    listagraficosRMSE[[len+1]] <- g1 + theme(legend.position = "none") 
    listagraficosACC[[len+1]] <- g2 + theme(legend.position = "none") 
    
    
  }# End Bin
}# End model

# Hago el grafico
# matrix 
legend <- LeyendaCondicional(Breaks = seq(-0.2,0.2,0.05), Paleta = "RdBu", Labels = c("MJO APORTA","MJO NO APORTA"))
t <- paste("pepe")
lg <- tableGrob(paste0(grupo,"-",model), theme= ttheme_minimal())
matrixgraficos <- rbind(matrix(1:(4*7),ncol = 4, nrow = 8, byrow = T),
                        c(8,8,8,8))
rg <- grid.arrange(grobs = listagraficos, ncol = 4, nrow = 8,
                   layout_matrix = matrixgraficos,
                   widths = rep(2.7,4),
                   heights = c(rep(2.5,7),0.2),
                   top = textGrob(t,gp=gpar(fontsize=13,font=3))) 

grid.arrange(bp, vp, legend, ncol=2, nrow = 2, 
             layout_matrix = rbind(c(1,2), c(3,3)),
             widths = c(2.7, 2.7), heights = c(2.5, 0.2))
grid.newpage()
grid.draw(cbind(lg, rg, size = "last"))

ggsave(file = OutFileName, arrangeGrob(grobs = plist, ncol = 2)) 
# Tabla de inicios activos e inactivos segun el modelo

