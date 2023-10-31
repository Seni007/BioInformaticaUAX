###############################################################################################
########################################## Trabajo R ##########################################
###############################################################################################

## Instalar RCurl
if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install("RCurl")


## Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/datos-trabajoR", followlocation = TRUE)
data = read.table(file = "datos-trabajoR.txt", head = TRUE)
##Para comprobar que la tabla es correcta
print(data)

###Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
###¿Cuántas variables hay? ¿Cuántos tratamientos? 
head(data)
summary(data)
dim(data)
str(data)
#Hay 2 variables. Hay 5 tratamientos 


###Haz un boxplot para nuestros datos. Uno para cada variable. 
###Colorea a Variable 1 y a Variable 2 de forma diferente (guarda esos colores para las siguientes gráficas)
boxplot(Variable1~Tratamiento, data=data, col = "#69e0fd", main = "BoxPlot Variable1")
boxplot(Variable2~Tratamiento, data=data, col = "#96fd69", main = "BoxPlot Variable2")



###Haz un gráfico de dispersión con las dos variables. Cada tratamiento debe de ir de un color distinto.
plot(x = data$Variable1, y = data$Variable2 , col = data$Tratamiento, 
     main = "Scatter plot de variables", xlab = "Variable1", ylab = "Variable2")



###Ponle leyenda al gráfico del apartado anterior. En el margen inferior derecho
legend("bottomright", legend = c("Tto1", "Tto2", "Tto3", "Tto4", "Tto5"), 
       fill = c("black", "darkred", "lightgreen", "turquoise", "blue"), title = "Legend")



###Haz un histograma para cada variable. Recuerda mantener los colores
hist(data$Variable1, col = "#69e0fd", xlab = "Variable 1", main = "Histogram Variable 1")
hist(data$Variable2, col = "#96fd69", xlab = "Variable 2", main = "Histogram Variable 2")



###Haz un factor en la columna tratamiento y guárdalo en una variable
factor_data_tto = factor(data$Tratamiento)



###Calcula la media y la desviación estándar para cada tratamiento.

##Usando aggregate, con cbind podemos promediar ambas columnas (Variable 1 y Variable 2)
resultado_aggregate = aggregate(cbind(Variable1, Variable2) ~ Tratamiento, data = data, mean)
##Renombro las columnas 
colnames(resultado_aggregate) = c("Tratamiento", "Media_Variable1", "Media_Variable2")
##Calcular la media de ambas variables y agregarla al resultado
resultado_aggregate$Media_Ambas = (resultado_aggregate$Media_Variable1 + resultado_aggregate$Media_Variable2) / 2
##Mostrar el resultado 
print(resultado_aggregate)

##Con tapply hacemos las medias de ambas variables 
media_variable1 = tapply(data$Variable1, data$Tratamiento, mean)
media_variable2 = tapply(data$Variable2, data$Tratamiento, mean)
##Hacemos la media de ambas variables
media_V1_V2 = (media_variable1 + media_variable2) / 2
##Crear un DataFrame con los resultados 
resultado_tapply <- data.frame(Tratamiento = unique(data$Tratamiento), 
                               Media_Variable1 = media_variable1, 
                               Media_Variable2 = media_variable2, 
                               Media_V1_V2 = media_V1_V2)
##Mostrar el resultado
print(resultado_tapply)
#Las medias de cada tratamiento son: Tto1: 2.255, Tto2: 

## Con aggregate para calcular la desviación estandar para cada tratamiento 
resultado_aggregate = aggregate(cbind(Variable1, Variable2) ~ Tratamiento, data = data, FUN = function(x) sd(x))
## Renombamos las columnas del resultado
colnames(resultado_aggregate) = c("Tratamiento", "SD_Variable1", "SD_Variable2")
## Calculamos la desviación estándar de ambas variables y agregarla al resultado
resultado_aggregate$SD_Ambas = sqrt(resultado_aggregate$SD_Variable1^2 + resultado_aggregate$SD_Variable2^2)
## Para ver el resultado
print(resultado_aggregate)

## Con tapply para calcular la desviación estandar para cada tratamiento.  
sd_variable1 = tapply(data$Variable1, data$Tratamiento, sd)
sd_variable2 = tapply(data$Variable2, data$Tratamiento, sd)
## Calculamos la desviación estándar de ambas variables
sd_ambas = sqrt(sd_variable1^2 + sd_variable2^2)
## Creamos un DataFrame con los resultados
resultado_tapply = data.frame(Tratamiento = unique(data$Tratamiento), 
                              SD_Variable1 = sd_variable1, 
                              SD_Variable2 = sd_variable2, 
                              SD_Ambas = sd_ambas)
## Para ver el resulado
print(resultado_tapply)


### Averigua cuántos elementos tiene cada tratamiento
numero_elementos_tratamiento = table(data$Tratamiento)
## Para mostrar la tabla que acabamos de hacer 
print(numero_elementos_tratamiento)


###Extrae los datos para el tratamiento 1 y el tratamiento 4 y guárdalos cada uno en una variable diferente.
bp = read.table("datos-trabajoR.txt",header=T)
Tratamiento_1 = data[data$Tratamiento == 1,]
Tratamiento_4 = data[data$Tratamiento == 4,]
print(Tratamiento_1)
print(Tratamiento_4)



### Nuestra hipótesis nula es que las medias de tratamiento 1 y tratamiento 4 para la Variable 1 son iguales.
### ¿Puedes comprobarlo? Para ello, necesitarás comprobar primero si los datos se distribuyen de forma normal.
### En función del resultado de la prueba de normalidad, ¿qué test usarías? 
###** En general, asumimos que las muestras son independientes, pero ¿son sus varianzas iguales? 
### Actúa de acuerdo a tus resultados.

## Usamos el shapiro test para ver si los datos se distribuyen de forma normal 
## Como vemos que si se distribuyen de forma normal (porque p>0.05) hacemos un t test para comprobar la hipotesis 
shapiro.test(Tratamiento_1$Variable1)
shapiro.test(Tratamiento_4$Variable1)

## Usamos el Fisher’s F-test para ver si nuestras dos muestras tienen la misma varianza 
## La hipotesis nula de es que las dos varianzas son iguales 
## Si el p value es mayor de 0.05 no podemos rechazar la hipotesis nula y las varianzas son iguales 
var.test(Tratamiento_1$Variable1, Tratamiento_4$Variable1)
# Como el p value es menor que 0.05 podemos rechazar la hipotesis nula y decir que las varianzas son DISTINTAS

## Ahora hacemos un t test para comprobar la hipotesis 
## Si el p value es mayor que 0.05 no podemos rechazar la hipotesis nula. 
## Lo que significa que los datos no son significativamente distintos 
t.test(Tratamiento_1$Variable1, Tratamiento_4$Variable1, var.equal = FALSE, paired = TRUE) # FALSE es si son distintos las varianzas
# Como el p value es menor que 0.05 podemos rechazar la hipotesis nula y los datos son significativamente distintos. 





