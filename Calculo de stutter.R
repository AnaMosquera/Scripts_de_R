rm(list = ls())
setwd("C:/Users/Ana/Desktop/Ana/Forense/2020_new STRs/Datos poblacionales/Stutter")

####################### CALCULO DEL PORCENTAJE DE STUTTER #############################

library(dplyr)

#Hay datos exportados del Genemapper en diferentes archivos, voy a leerlos y unirlos:
#Hay que SACAR estos archivos con los datos del DYE

Datos1 <- read.table("1STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos1)
Datos1 <- Datos1[,-7]#Aparece una columna que no tiene nada, la quito.
head(Datos1)

Datos2 <- read.table("2STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos2)
Datos2 <- Datos2[,-7]
head(Datos2)

Datos3 <- read.table("3STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos3)
Datos3 <- Datos3[,-7]
head(Datos3)

Datos4 <- read.table("4STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos4)
Datos4 <- Datos4[,-7]
head(Datos4)

Datos5 <- read.table("5STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos5)
Datos5 <- Datos5[,-7]
head(Datos5)

Datos6 <- read.table("6STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos6)
Datos6 <- Datos6[,-7]
head(Datos6)

Datos7 <- read.table("7STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos7)
Datos7 <- Datos7[,-7]
head(Datos7)

Datos8 <- read.table("8STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos8)
Datos8 <- Datos8[,-7]
head(Datos8)

Datos9 <- read.table("9STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos9)
Datos9 <- Datos9[,-7]
head(Datos9)

Datos10 <- read.table("10STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos10)
Datos10 <- Datos10[,-7]
head(Datos10)

Datos11 <- read.table("11STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos11)
Datos11 <- Datos11[,-7]
head(Datos11)

Datos12 <- read.table("12STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos12)
Datos12 <- Datos12[,-7]
head(Datos12)

Datos13 <- read.table("13STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos13)
Datos13 <- Datos13[,-7]
head(Datos13)

Datos14 <- read.table("14STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos14)
Datos14 <- Datos14[,-7]
head(Datos14)

Datos15 <- read.table("15STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos15)
Datos15 <- Datos15[,-7]
head(Datos15)

Datos16 <- read.table("16STRs nuevos_%stutter_SamplePlotSizingTable.txt",header = T, sep = "\t")
head(Datos16)
Datos16 <- Datos16[,-7]
head(Datos16)

Datos<-rbind(Datos1,Datos2,Datos3,Datos4,Datos5,Datos6,Datos7,Datos8,Datos9,Datos10,Datos11,Datos12,Datos13,Datos14,Datos15,Datos16)
summary(Datos)

Datos$Marker <- factor(Datos$Marker)
Datos$Allele <- factor(Datos$Allele)
levels(Datos$Marker)
summary(Datos)
head(Datos)

#Para completar las celdas vacías de Marker como NA:
is.na(Datos$Marker)#no los detecta como NA
Datos$Marker <- ifelse(Datos$Marker == "", NA, Datos$Marker)#sustituto el level "" por NA, y lo demás queda igual
is.na(Datos$Marker)
levels(Datos$Marker)
Datos$Marker <- factor(Datos$Marker)
levels(Datos$Marker)#No sé pq borra los nombres y les pone números, así que se los vuelvo a asignar
cromosomas<-c("Chr1","Chr13","Chr15","Chr19","Chr20","Chr3")
levels(Datos$Marker)<-cromosomas
levels(Datos$Marker)
summary(Datos)#Son los mismos que antes, así que está bien
head(Datos)

table(Datos$Marker,Datos$Allele)
write.table(Datos,file="Datos_Stutter.txt", col.names=T,row.names=F)

#En la tabla "Datos_Stutter" resultante, en la columna "Marker", las celdas vacías serían picos stutter que pueden tener una unidad de repetición menos o más que los alelos detectados.
#El objetivo es calcular el % de stutter plus y minus de cada alelo, definiendo el stutter plus como aquellos alelos con un tamaño 5pb mayor que el alelo y el stutter minus aquellos alelos con un tamaño de 5pb menor al alelo, SOLO cuando las celdas de "Marker" y "Allele" estén vacías.

#Creo una variable nueva que sea el valor de Size redondeado (para que esas diferencias de 5pb se basen en SizeRound):
Datos <- Datos %>% mutate (SizeRound=round(Size,digits = 0),.after = Size)

#Ordeno por Size, para cada muestra:
Datos <- Datos %>% arrange(Sample.File.Name,Size)

#Creo dos nuevas columnas, con NAs:
Datos <- Datos %>% mutate(StutterMinus=NA,StutterPlus=NA)

#Para cada muestra, quiero que, si para la fila [i] el Marker=NA y Allele=NA, mire la fila siguiente [i+1] y si SizeRound [i+1]-[i] es 5pb, calcule el % de stutter, copie el Marker de la fila [i+1] y añada en Allele "Stutter Minus":

Datos_completados <- Datos #Creo una copia del dataframe original
#La fila 2, que me parece rara, no sé si será fluorescencia residual (REVISAR), la saco:
Datos_completados <-Datos_completados[-2,]

#"lag" hace referencia a la fila anterior y "lead" a la siguiente

Datos_completados <- Datos_completados %>% mutate(StutterMinus=ifelse(
  Sample.File.Name==lag(Sample.File.Name) & is.na(lag(Marker)) & !is.na(Marker) & is.na(lag(Allele)) & !is.na(Allele) & Height>3*lag(Height) & (SizeRound-lag(SizeRound)==5 |round(Size-lag(Size)==5)), lag(Height)/Height*100,NA))
#Así añado la condición de que SizeRound-lag(SizeRound)==5 O round(Size-lag(Size)==5)

summary(Datos_completados$StutterMinus)
MediaStutterMinus<-mean(Datos_completados$StutterMinus,na.rm = T)
MediaStutterMinus

indice_max_StutterMinus <- which.max(Datos_completados$StutterMinus)
indice_max_StutterMinus
cat("El valor máximo en la columna está en la línea:", indice_max_StutterMinus, "\n")

#Había un valor máximo de 2192, pq no había puesto de condición que la fila del alelo no tenía NA en Marker y Allele. Corregido.

#Luego había otro max de 320, en la fila 132, pq coincidía que 2 alelos, de diferentes marcadores, tienen un tamaño muy similar y parece que es stutter de la fila [i+2](Chr19), no de la [i+1] (Chr1)... Añado la condición de que el la altura de stutter tiene que ser al menos 3 veces menor de la del alelo. Así no calcula con el siguiente marcador, pero al menos no da un % incorrecto.

#El siguiente max es de 26, en las líneas 245-246. Si nos fijamos en los decimales de Size, casi serían 6pb... Además, hay 2 líneas en las que el Size es exactamente el mismo (289.72), pero con valores de altura diferentes. Creo que no es su stutter. (REVISAR)

Top10_StutterMinus <- Datos_completados %>% arrange(desc(StutterMinus)) %>% head(10) # Ordena el dataframe de forma descendente por la columna de interés, y obtiene los 10 valores primeros (para todas las columnas)
Top10_StutterMinus
#De estos valores, hay que revisar por lo menos el de 20% y 18% pq podría ocurrir lo mismo

valor_buscar <- "C10_C10_PLACA178.hid"  # Reemplazar el nº a buscar
filas_con_valor <- which(Datos_completados$Sample.File.Name == valor_buscar)
if (length(filas_con_valor) > 0) {
  cat("El valor", valor_buscar, "se encuentra en las siguientes líneas:", filas_con_valor, "\n")
} else {
  cat("El valor", valor_buscar, "no se encontró en ninguna línea.\n")
}#Saca una frase hecha en la pantalla

#No sé pq (igual por diferencias en decimales) no consigo localizar el valor de 20% de stutter minus usando esa columna, así que lo busco por muestra (linea 1375). Pasa lo mismo que con el valor de 26%, hay dos filas con el mismo valor 299.79.

#Lo mismo con 18% (linea 82-83), dos filas con Size 289.78

#Para los valores de 15% y 14% de stutter minus ya no ocurre eso, COMPROBAR que no es fluorescencia residual pero puede que estos ya sean valores reales de stutter

#Creo una columna extra: "Revisar" y pongo el motivo: 

Datos_completados <- Datos_completados %>% mutate (Revisar=NA)

Datos_completados <-Datos_completados %>% mutate(Revisar=ifelse(Sample.File.Name==lag(Sample.File.Name) & is.na(Marker) & !is.na(lag(Marker)) & is.na(Allele) & !is.na(lag(Allele)) & Size==lag(Size) & lag(Height)>8*Height,"Fluorescencia?",NA))#Puede ser fluorescencia residual de otro color

Fluorescencia <- which(Datos_completados$Revisar == "Fluorescencia?")
Fluorescencia #no sé pq en la visualización los veo como un nº más que el que reporta como la fila donde aparece "Fluorescencia?"

Datos_completados <- Datos_completados %>% mutate(Revisar=ifelse(Sample.File.Name==lag(Sample.File.Name) & !is.na(StutterMinus) & lag(Revisar=="Fluorescencia?"), "StutterM_incorrecto",Revisar))#Indico qué valores de stutter no me creo, por fluorescencia

StutterMinus_incorrecto <- which(Datos_completados$Revisar=="StutterM_incorrecto")
StutterMinus_incorrecto #Yo veo como una fila más, pero no es real 

#Elimino los % de stutter que no nos creemos
filas_a_modificar <- c(83, 245, 989, 1109, 1374, 1905, 1914, 2153, 2332, 2707)
Datos_completados[filas_a_modificar, "StutterMinus"]<-NA

#Si se confirma que es fluorescencia residual, se podrían eliminar esas filas (6,82,229,244,254,303,330,385,416,428,555,988,1108,1176,1373,1745,1779,1904,1913,1986,2003,2152,2331,2706,2824,2852) y así probablemente podríamos calcular el stutter minus de algún dato más.

summary(Datos_completados$StutterMinus)
MediaStutterMinus<-mean(Datos_completados$StutterMinus,na.rm = T)
MediaStutterMinus


#Para hacer lo mismo con el Stutter Plus:
Datos_completados <- Datos_completados %>% mutate(StutterPlus=ifelse(
  Sample.File.Name==lead(Sample.File.Name) & is.na(lead(Marker)) & !is.na(Marker) & is.na(lead(Allele)) & !is.na(Allele) & Height>3*lead(Height) & (lead(SizeRound)-SizeRound==5 |round(lead(Size)-Size)==5), 
  lead(Height)/Height*100,NA))

summary(Datos_completados$StutterPlus)
MediaStutterPlus<-mean(Datos_completados$StutterPlus,na.rm = T)
MediaStutterPlus

indice_max_StutterMinus<- which.max(Datos_completados$StutterPlus)
indice_max_StutterMinus
#Hay un valor máximo de 32, y tiene que ver con que calcula el stutter plus de un stutter minus, creo que de diferente dye,... hay que incluir los dyes en el los datos de partida (REVISAR)...

Top10_StutterPlus <- Datos_completados %>% arrange(desc(StutterPlus)) %>% head(10)
Top10_StutterPlus

texto_a_bucar <- "L02_placa2_6STRs_nuevos_G08_PLACA177.hid"
filas_valor <- which(Datos_completados$Sample.File.Name == texto_a_bucar)
filas_valor
#REVISAR LÍNEAS: 40

#OJO: stutter Plus de un alelo se suma al stutter minus de otro alelo (mismo Chr)
Datos_completados <- Datos_completados %>% mutate(Revisar=ifelse(
  Sample.File.Name==lag(Sample.File.Name) & Sample.File.Name==lead(Sample.File.Name) & 
    is.na(Marker) & is.na(Allele) & 
    !is.na(lag(Marker)) & !is.na(lag(Allele)) & !is.na(lag(StutterPlus)) &
    !is.na(lead(Marker)) & !is.na(lead(Allele)) & !is.na(lead(StutterMinus)) & 
    lag(Marker)==lead(Marker) & lead(SizeRound)-lag(SizeRound)==10
    ,"Chr_stutter_M+P",Revisar))

Stutter_M_P_combinado <-which(Datos_completados$Revisar=="Chr_stutter_M+P")
Stutter_M_P_combinado

summary(Datos_completados$StutterMinus)
summary(Datos_completados$StutterPlus)

filas_a_modificar_StutterPlus <- c(17, 100, 177,  188,  312,  451,  459,  471,  562,  989, 1053, 1099, 1138, 1367, 1374, 1721, 1967, 2067, 2232, 2274, 2428, 2579, 2629)
Datos_completados[filas_a_modificar_StutterPlus, "StutterPlus"]<-NA

filas_a_modificar_StutterMinus <- c(19, 102, 179,  190,  314,  453,  461,  473,  564,  991, 1055, 1101, 1140, 1369, 1376, 1723, 1969, 2069, 2234, 2276, 2430, 2581, 2631)
Datos_completados[filas_a_modificar_StutterMinus, "StutterMinus"]<-NA

summary(Datos_completados$StutterMinus)
summary(Datos_completados$StutterPlus)

write.table(Datos_completados,file="Datos_Stutter_completados.txt", col.names=T,row.names=F)

#Calcula las medias de Stutter Minus y Plus por cada nivel de Marker
StutterMinus_por_marcador <- tapply(Datos_completados$StutterMinus,Datos_completados$Marker,FUN = mean,na.rm = T)
StutterMinus_por_marcador

StutterMinus_SD_por_marcador <- tapply(Datos_completados$StutterMinus,Datos_completados$Marker,FUN = sd,na.rm = T)
StutterMinus_SD_por_marcador

StutterMinus_Mean_3SD_por_marcador <- StutterMinus_por_marcador+3*StutterMinus_SD_por_marcador
StutterMinus_Mean_3SD_por_marcador


StutterPlus_por_marcador <- tapply(Datos_completados$StutterPlus,Datos_completados$Marker,FUN = mean, na.rm=T)
StutterPlus_por_marcador

StutterPlus_SD_por_marcador <- tapply(Datos_completados$StutterPlus,Datos_completados$Marker,FUN = sd,na.rm = T)
StutterPlus_SD_por_marcador 

StutterPlus_Mean_3SD_por_marcador <- StutterPlus_por_marcador+3*StutterPlus_SD_por_marcador
StutterPlus_Mean_3SD_por_marcador


#PENDIENTE:

#Correr el script (y adaptarlo) con los datos de los Dyes

#Revisar datos de fluorescencia residual, y eliminarlos si se confirma

#Revisar cuando el stutter Plus de un alelo se suma al stutter minus de otro alelo (diferente Chr, diferente Dye)

# ¿Rehacer filtrando por el umbral analítico? Dado que hemos establecido el umbral analítico (UA) en 180pb
# HACER un nuevo script con este filtro y adaptarlo
#min(Datos$Height)
#DatosUA <- subset(Datos,Height>=180)
#min(DatosUA$Height)
