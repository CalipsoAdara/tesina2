# Analisis de la ola de calor 2013. 11/12 al 31/12
#
# Basado en el paper de los dioses
#Subseasonal prediction of the heat wave of December 2013
#in Southern South America by the POAMA and BCC‐CPS models
#
# By lucia m. castro
# ---------------------------------------------------------------------------------

# Limpiar el espacio
rm(list=ls())

# cargar paquetes
library(ggplot2)
library(metR)

# Cargo mis funciones
source("/vegeta/datos/SubX/mjo.lucia.cas/funciones.R")
svpath = "/home/lucia.castro/tesina2/resultados"
setwd("/vegeta/datos/SubX/mjo.lucia.cas")

# Modelos
groups=c('GMAO','RSMAS','ESRL','ECCC','NRL','EMC','MME')                       
models=c('GEOS_V2p1','CCSM4','FIMr1p1','GEM','NESM','GEFS','SAT') 
glabel=paste0(groups,"-",models)

# Dias de la ola de calor
diasola = seq.Date(as.Date("2013-12-13"),as.Date("2013-12-31"),1)
olaw1 = diasola[1:7]  # 13 a 19
olaw2 = diasola[8:14] # 20 a 26 
olaw3 = diasola[15:19]  # 5 dias
olaw <- list(olaw1,olaw2,olaw3)
leadbin <- c(7,14,21)


# Primera forma: Analizar como los modelos pronosticaro cada semana de la ola. 
# todos los modelos tienen inicializaciones en fechas diferentes, por lo tanto el lead a las
# semanas de la ola seran distintos para c/u. VOy a considerar grupos de leads por ej (leads 2-5,lead 5-8)
listagraficos <- list()

for(w in olaw) { # por cada semana del evento
  for(m in groups) { # por cada modelo
    
    # Cargar datos del modelo (solo el subset nov-dic de 2013)
    targetdate <- readRDS(paste0("targetdate_",m,"_ONDEFM.rds"))
    startdate <- dimnames(targetdate)$startdate
    mod <- readRDS(paste0("model_olacalor_",m,".rds"))
    
    for(l in leadbin) { # por cada lead 
      # Buscar la inicializacion anterior a la semana del evento
      sem_anterior = as.character(w-l)
      
      # CONDICIONAL: en la semana 3 solo hay 5 dias, extiendo "sem_anterior" para tener en cuenta los 7 lead posibles
      if(length(w) == 5 ) {
        date_lead_mayor = (w-l)[1]
        date_lead_menor = (w-l)[1] +6
        sem_anterior = as.character(seq.Date(date_lead_mayor,date_lead_menor,by = 1))}
      
      stdate_ola <- last(intersect(startdate, sem_anterior))
      
      # Osea que el lead al evento es:
      lead =  as.numeric(w[1] - as.Date(stdate_ola))
      
      # Buscar las fechas de la semana hindcasteadas 
      hind_ola <- targetdate[,stdate_ola] %in% as.character(w)
      
      # Evaluo en el stdate de la ola y los dias de la semana. luego promedio
      mod_ola <- apply(mod[,,hind_ola,stdate_ola], c(1,2), FUN = mean)
      g <- GraphDiscrete(Data = ggDataFrame(mod_ola),Breaks = seq(-3,3,1),Label = "°C",
                    Paleta = "RdBu", Direccion = -1, 
                    Titulo = paste0("LEAD ", l-6,"-",l)) + labs(title = NULL)
      
      # Lo completo asi en vez de append porque sino guarda listas en vez de ggplots
      len <- length(listagraficos)
      listagraficos[[len+1]] <- g  
      
    }# End loop lead
  } # End loop modelo
} # End loop semana


# G R A F I C O S ---------------------------------------------------------
# uno por cada semana

# De esta forma puedo poner los nombres de los modelos en las filas
png(paste0(svpath,"/olaw2.png"),width = 24, height = 30, units = "cm",res = 1080) # creo el objeto

# lg <- tableGrob(c("", glabel), theme= ttheme_minimal(base_size = 8))     # Nombre fila
# rg <- arrangeGrob(grobs = listagraficos[1:(7*3)], ncol=3,
#                   top = textGrob("Ola de calor Dic 2013 \n 13/12 - 19/12"
#                                  ,gp=gpar(fontsize=18)))


fil = c("LEAD 1-7","LEAD 8-14","LEAD 15-21")
combine <- rbind(tableGrob(t(fil), theme = ttheme_minimal(), rows = ""), 
                 cbind(tableGrob(glabel, theme = ttheme_minimal()), 
                       arrangeGrob(grobs = listagraficos[22:42],nrow = 7,ncol=3),  size = "last"), size = "last")

grid.newpage()
grid.draw(combine)

dev.off() # borra


# Segunda forma: Analizar las medias semanales ya calculadas de los mod vs las obs
# 

listagraficos2 <- list()
for(mod in groups) { # por cada modelo
  #cargo el targatdate y el startdate
  targetdate <- readRDS(paste0("targetdate_",mod,"_ONDEFM.rds"))
  startdate <- dimnames(targetdate)$startdate
  
  # Busco que dias de la ola fueron modeladas. 
  # Tomo la el startdate que caiga en la semana anterior a la ola (si hay dos eligo el ultimo)
  stdate_ola <- last(intersect(startdate, sem_anterior))
  
  # Cargo los datos semanales de los modelos
  modweek <- readRDS(paste0("modelweek_",mod,".rds"))
  obsweek <- readRDS(paste0("obsweek_",mod,".rds"))
  
  # Evaluo en la fecha y resto
  resta = obsweek[,,stdate_ola,] - modweek[,,stdate_ola,]
  # Grafico
  g <- GraphDiscreteMultiple(Data = ggScoreSemanal(resta),Breaks = seq(-5,5,1),
                             Label = "°C",Paleta = "RdBu", Direccion = -1)
  # Lo completo asi en vez de append porque sino guarda listas en vez de ggplots
  len <- length(listagraficos2)
  listagraficos2[[len+1]] <- g  
  
} # End model loop  

# De esta forma puedo poner los nombres de los modelos en las filas
png(paste0(svpath,"/ola.png"),width = 20, height = 35, units = "cm",res = 1080) # creo el objeto

lg <- tableGrob(c("", glabel), theme= ttheme_minimal(base_size = 10))
rg <- arrangeGrob(grobs = listagraficos2, ncol=1,
                  top = textGrob("Ola de calor Dic 2013 \n OBS-MOD"
                                 ,gp=gpar(fontsize=18)))

grid.newpage()
grid.draw(cbind(lg, rg, size = "last")) # lo dibuja

dev.off() # borra

# ------------------------------------------------------------------------------

# O B S E R V A C I O N E S 

# Cargo datos
ar.anom = readRDS("/datos4/Obs/t2m_cpc_daily/t2manom_NOAA.rds")
caso = 0
for(w in olaw) {
  
  sem <- format(w, "%d/%m")
  # Promedio semanal
  prom <- apply(ar.anom[,,as.character(w)],c(1,2),FUN = mean)
  # Graficar
  go <- GraphDiscrete(Data = ggDataFrame(prom),Breaks = seq(-3,3,1),Label = "°C",
                           Paleta = "RdBu", Direccion = -1, 
                           Titulo = paste0("T2M ANOM CPC \n",sem[1],"-",sem[length(sem)]))
  caso <- caso +1 
  ggsave(paste0(svpath,"/obsola_",caso,".png"),go)
  
}

