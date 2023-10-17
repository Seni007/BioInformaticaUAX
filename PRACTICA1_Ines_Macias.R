#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 03 OCTUBRE 23:59
## Se requiere la entrega de este script completado con los códigos más las imágenes y las respuestas a las preguntas
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) 
head(data)
tail(data)


# Hacemos un primer histograma para explorar los datos
hist(data)
hist(data, col = "gray", main="GSE5583 - Histogram")


# Transformamos los datos con un logaritmo 
# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
data_log=log2 (data) #Transformamos los datos en otro variable. Se guarda con el nombre de data_log para sacar un histograma donde se vean los datos mejor con logaritmo base 2. Ahora los datos se ven con forma de campana de gaus.
hist(data_log) #Hacemos el histograma de el logaritmo de los datos


# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
# ¿Qué es un boxplot?
boxplot(data_log) #Box plot sin colores 
boxplot (data_log, col=c("blue","blue","blue","orange","orange","orange"), main="GSE5583 - boxplots", las=2) #Los colores se ponen con col=c, los wild types estan en azul y los knock out en naranja. El titulo se pone con main=. Las=2 es para tener los ejes en vertical.
#Un diagrama de cajas es una técnica convencional utilizada para visualizar una serie de valores numéricos mediante la representación gráfica de sus cuartiles. De este modo, se hace evidente de manera sencilla la mediana y los valores que delimitan los cuartiles en los datos, y también se pueden identificar los datos atípicos.


# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
# de los valores de expresión. ¿Es correcta la separación?
hc = hclust (as.dist(1-cor(data_log))) 
plot(hc, main="Hierarchical Clustering") #La separación es correcta (hay 3 KO y 3 WT)
#Hierarchical clustering es un metodo de analisis en el que agrupaciones por tipos 


#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
wt <- data[,1:3] #WT van de la columna 1-3. La coma es porque solo se quiere que se cojan las columnas.
ko <- data[,4:6] #KO van de columna 4-6. 
class(wt)
#Hemos generado un matrix array, una tabla 


# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean = apply(wt, 1, mean) #Se calcula la media de todas las filas con la función apply. Creamos una variable para tener las medias, en este caso wt.mean. El 1 es para calcular la media para cada fila (2 sería para cada columna).
head(wt.mean) 
ko.mean = apply(ko, 1, mean) #Lo mismo se hace con los KO
head(ko.mean)


# ¿Cuál es la media más alta?
max (wt.mean)
max (ko.mean)
#La más alta es el KO


# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean ~ wt.mean) #Estamos enfrentando las dos muestras
plot(ko.mean ~ wt.mean, xlab = "WT", ylab = "KO", main = "GSE5583 - Scatter") #xlab es para cambiar el lable del x-axis, ylab es para cambiar el lable del y-axis. Main para ponerle el titulo.


# Añadir una línea diagonal con abline
abline(0, 1, col = "red") #b=0 y=1*x
abline(h=2, col = "blue") #La "h" es para que nos de una linea horizontal 
abline(v=5, col = "green") #"v" es para ina linea vertical 

# ¿Eres capaz de añadirle un grid?
grid()


# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean
diff.mean


# Hacemos un histograma de las diferencias de medias
hist(diff.mean) #Sin logaritmo sale un poco mal 
diff.mean_log=log2 (diff.mean) #guardamos el logaritmo de diff.mean como diff.mean_log
hist(diff.mean_log, col = "turquoise")



# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.
# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué?
# ¿Cuántas valores tiene cada muestra?
pvalue = NULL
tstat = NULL
for(i in 1 : nrow(data)) { #para cada gen
  x = wt[i,] # gene wt número i
  y = ko[i,] # gene ko número i
  
  # Hacemos el test 
  t = t.test(x, y)
  
  # Añadimos el p-value a la lista
  pvalue[i] = t$p.value
  # Añadimos las estadísticas a la lista 
  tstat[i] = t$statistic
  }
# "i" es solo una variable que existe en el bucle 
# wt y ko son vectores. x,y son vectores. El t test se guarda en una variable que llamamos "t"

head(pvalue)
length(pvalue)

# Ahora comprobamos que hemos hecho TODOS los cálculos


# Hacemos un histograma de los p-values.
# ¿Qué pasa si le ponemos con una transformación de -log10?
hist(pvalue)
hist(-log10(pvalue), col = "purple")


# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano") #Tenemos la diferencia de media contra el log-10 del p value
# Los p values significativos son los que están por arriba del volcano 


# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff = 2 #la diferencia de medias está en 2 por arriba y -2 abajo
pvalue_cutoff = 0.01
abline(v = diff.mean_cutoff, col = "blue", lwd = 3) # Abline nos da una recta, lwd es la anchura de la linea 
#abline(v = -diff.mean_cutoff, col = "blue", lwd = 3) este seria el -2 pero va a solapar con el +2
abline(h = -log10(pvalue_cutoff), col = "green", lwd = 3) # Hay que convertirlo a log porque sin o no el codigo no lo va a reconocer 


# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean = abs(diff.mean) >= diff.mean_cutoff
dim(data[filter_by_diff.mean, ])
# Filtro para la diferencia de medias, 2 y -2. Abs es para dar un valor absoluto, sacamos todos los que seam por encima de 2. 


# Ahora el filtro de p-value
filter_by_pvalue = pvalue <= pvalue_cutoff #Extraemos tosos los pvalues que sean iguales o menores que el p value que hemos definido 
dim(data[filter_by_pvalue, ]) # Dim es para las dimensiones 


# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios?
filter_combined = filter_by_diff.mean & filter_by_pvalue
filtered = data [filter_combined,] # Combinamos ambos identificadores y nos quedamos solo con los genes que coincidan en los dos filtros 
dim(filtered) 
head(filtered)


# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #2")
points (diff.mean[filter_combined], -log10(pvalue[filter_combined]), col = "red") 


# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés?
plot(diff.mean, -log10(pvalue), main = "GSE5583 - Volcano #3")
points (diff.mean[filter_combined & diff.mean < 0],
        -log10(pvalue[filter_combined & diff.mean < 0]), col = "red")
points (diff.mean[filter_combined & diff.mean > 0],
        -log10(pvalue[filter_combined & diff.mean > 0]), col = "blue")
#diff.mean= wt.mean -ko.mean


# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
heatmap(filtered) # Sin ordenar

rowv = as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv = as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv,Colv=colv, cexCol=0.7,labRow=FALSE) # Los hemos agrupado, los wt por un lado y los ko por otro 


# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)
library(RColorBrewer)

# Hacemos nuestro heatmap
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
          col = rev(redblue(256)), scale = "row")


# Lo guardamos en un archivo PDF
pdf ("GSE5583_DE_Heatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
          col= rev(redblue(256)), scale = "row",labRow=FALSE)
dev.off()
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
          col = rev(redgreen(75)), scale = "row",labRow=FALSE)



# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table (filtered, "GSE5583_DE.txt", sep ="\t",
              quote = FALSE)

