if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")
if (!requireNamespace("limma", quietly = TRUE))
  install.packages("limma")
if (!requireNamespace("SummarizedExperiment", quietly = TRUE))
  install.packages("SummarizedExperiment")
if (!requireNamespace("readxl", quietly = TRUE))
  install.packages("readxl")

library(SummarizedExperiment)
library (ggplot2)
library (dplyr)
library (readxl)
library (limma)

datos <- read_excel("C:/Users/Victor/Desktop/R/datos/TIO2+PTYR-human-MSS+MSIvsPD.xlsx")
head (datos)

matriz_datos <- as.matrix(datos[,c("M1_1_MSS","M1_2_MSS","M5_1_MSS","M5_2_MSS", "T49_1_MSS", "T49_2_MSS", "M42_1_PD", "M42_2_PD", "M43_1_PD", "M43_2_PD", "M64_1_PD", "M64_2_PD")])

Metadata <- data.frame (sample_id = colnames(matriz_datos), condition = c(rep("MSS", 6), rep ("PD", 6)), strigsAsFactors = FALSE)

CaracteristicasMetadata <- data.frame(
  SequenceModifications = datos$SequenceModifications, 
  stringsAsFactors = FALSE
)
CaracteristicasMetadata

#Creamos el objeto Summarized Experiment
se <- SummarizedExperiment(
  assays = list(counts = matriz_datos),
  rowData = CaracteristicasMetadata,
  colData = Metadata,
)
#Comprobamos el Summarized Experiment
assay(se, "counts")
colData(se)
rowData(se)
head (se$)

#Obtenemos las dimensiones dle dataset
dim(se)

#Representamos la distribución de los datos
boxplot (assay(se, "counts"), main = "Distribución de los datos", xlab = "muestras", ylab = "niveles de expresión", las = 2, col = "lightblue")


#Representamos el histograma de las 5 primeras muestras para ver su distribución más detalladamente 
par (mfrow = c(1, 3))
for (i in 1:5) {
  hist(assay(se, "counts")[, i], main = paste("Histograma de", colnames(se)[i]),
       xlab = "Abundancia", col = "lightgreen", breaks = 30)
}

#Transformamos los datos en base de logaritmo base 2, sumamos uno ya que algunos valores son 0 y de esta forma evitamos errores en la transformación
log_se <- log2(assay(se, "counts") + 1)
#Representamos los datos para ver la nueva distribución de los datos
colores <- c(rep("red", 2), rep("blue", 2), rep("green", 2), rep("orange", 2), rep("grey", 2), rep("lightblue", 2))
boxplot (log_se, main = "Distribución logaritmica", xlab = "Muestras", ylab = "log2(Expresión", las = 2, col = colores)
par (mfrow = c(1, 5))
for (i in 1:5) {
  hist(log_se[, i], main = paste("Histograma de", colnames(log_se)[i]),
       xlab = "Abundancia", col = "lightgreen", breaks = 30)
}

#Calculamos la matriz de correlaciones
Cor_se <- cor (log_se)

#Y la visualizamos mediante un mapa de distancias 
heatmap (Cor_se, sym = TRUE, col = colorRampPalette(c("blue", "white", "red"))(100), main = "Matriz de correlación")
legend ("topright",
        legend = c("Correlación negativa", "Sin correlación", "Correlación positiva"),
        fill = c("blue", "white", "red"),
        title = "correlación",
        cex = 0.5)

PCA <- prcomp(t(log_se), scale = FALSE)
loads <- round (PCA$sdev^2/sum (PCA$sdev^2)*100, 1)
#Visualizamos el PCA
xlab <- c(paste("PCA1", loads[1], "%"))
ylab <- c(paste ("PCA2", loads [2], "%"))
plot(PCA$x[, 1:2], col = colores, pch = 19, xlab = xlab, ylab = ylab, main ="PCA")

legend ("topright", 
        legend = unique(Metadata$sample_id),
        col = colores,
        pch = 19,
        title = "Condiciones",
        cex = 0.6)
#Con este gráfico no podemos obtener demasiada información, más allá de que las condiciones M42_1_PD y M42_2_PD están bastante aleadas del resto en el gráfico. Representamos un gráfico tridimensional con PCA3 para ver si podemos obtener más información

library(scatterplot3d)

scatterplot3d(
  PCA$x[, 1],
  PCA$x[, 2],
  PCA$x[, 3],
  color = colores, 
  pch = 19,
  xlab = paste("PCA1", loads[1], "%"),
  ylab = paste("PCA2", loads[2], "%"),
  zlab = paste("PCA3", loads[3], "%"),
  main = "PCA en 3D"
)
legend ("topright",
        legend = unique(Metadata$sample_id),
        col = unique(colores),
        pch = 15,
        title = "Condiciones",
        cex = 0.6)

#T test
t_test <- apply (log_se, 1, function(x){
  grupo_MSS <- x[Metadata$condition == "MSS"]
  grupo_PD <- x[Metadata$condition == "PD"]
  t.test(grupo_MSS, grupo_PD)$p.value
})

#calculamos cuantos de los valores obtenidos son significativos
suma <- sum ((t_test) < 0.05, na.rm = TRUE)
suma

#Creamos un data frame con los resultados significativos
p_value_Significativo <- data.frame (
  Fosfopeptido = CaracteristicasMetadata$SequenceModifications,
  p_value = t_test < 0.05,
  valores = t_test
)

# Calcular el log2 del Fold Change
logFC <- rowMeans(log_se[, Metadata$condition == "PD"], na.rm = TRUE) - rowMeans(log_se[, Metadata$condition == "MSS"], na.rm = TRUE)




# Crear un data frame para el Volcano Plot
volcano_data <- data.frame(
  Fosfopeptido = CaracteristicasMetadata$SequenceModifications,
  id = datos$Accession,
  logFC = logFC,
  p_value = t_test,
  negLogP = -log10(t_test)
)

# Agregar la columna de significatividad
volcano_data$Significativo <- ifelse(volcano_data$p_value < 0.05, "Significativo", "No Significativo")

# Filtrar solo los fosfopeptidos significativos
volcano_data_significativos <- volcano_data[volcano_data$Significativo == "Significativo", ]

# Generar el Volcano Plot
library(ggplot2)

ggplot(volcano_data, aes(x = logFC, y = negLogP, color = Significativo)) + 
  geom_point(alpha = 0.7) + 
  scale_color_manual(values = c("No Significativo" = "grey", "Significativo" = "red")) + 
  labs(title = "Volcano Plot de Fosfopeptidos Significativos",
       x = "Log2 Fold Change",
       y = "-Log10(p-valor)") + 
  theme_minimal()

# Añadir información sobre el tipo de cáncer
volcano_data$Tipo_Cancer <- ifelse(volcano_data$logFC > 0, "PD", "MSS")

# Filtrar solo los fosfopeptidos significativos
volcano_data_significativos <- volcano_data[volcano_data$Significativo == "Significativo", ]

# Mostrar los resultados con el tipo de cáncer característico
print(volcano_data_significativos[, c("Fosfopeptido", "logFC", "p_value", "Tipo_Cancer", "id")])

# Diseñamos el modelo lineal
diseño <- model.matrix(~ condition, data = Metadata)

# Ajustamos el modelo
fit <- lmFit(log_se, diseño)
fit <- eBayes(fit)

# Obtenemos los resultados de los fosfopéptidos
resultados <- topTable(fit, coef = "conditionPD", adjust = "BH", number = Inf)
resultados$Fosfopeptido <- CaracteristicasMetadata$SequenceModifications
resultados$id <- datos$Accession


#Data frame
data_resultados <- data.frame(
  Fosfopeptido = CaracteristicasMetadata$SequenceModifications,
  id = datos$Accession,
  logFC = resultados$logFC,
  p_value_adj = resultados$adj.P.Val,
  negLogP = -log10(t_test)
)

resultados$Significativo <- ifelse(resultados$adj.P.Val < 0.05 & abs(resultados$logFC) > 1, "Significativo", "No Significativo")

# Añadir información sobre el tipo de cáncer a los fosfopeptidos diferenciales
data_resultados$Tipo_Cancer <- ifelse(resultados$logFC > 0, "PD", "MSS")

# Filtramos los fosfopéptidos con valor p ajustado significativo
fosfopeptidos_diferenciales <- data_resultados %>% filter(p_value_adj < 0.05)

# Contamos los fosfopéptidos significativos
num_significativos <- nrow(fosfopeptidos_diferenciales)
cat("Número de fosfopéptidos significativos:", num_significativos, "\n")

# Preparación del volcano plot
resultados$logFC <- resultados$logFC
resultados$negLogP <- -log10(resultados$adj.P.Val)

print (fosfopeptidos_diferenciales)

#Generamos el Volcano plot 
ggplot (resultados, aes (x=logFC, y = negLogP, color = Significativo)) + geom_point(alpha = 0.7) + scale_color_manual(values = c("No significativo" = "grey", "Significativo" = "red")) + labs(title = "Volcano Plot de fosfopeptidos significativos (logaritmico)", x = "Log2", y = "Log2") + theme_minimal()


# Obtener las modificaciones de secuencia
mod_fosfopeptidos <- data.frame(Fosfopeptido = fosfopeptidos_diferenciales$Fosfopeptido, 
                                Tipo_Cancer = fosfopeptidos_diferenciales$Tipo_Cancer)
mod_volcano_plot <- data.frame(Fosfopeptido = volcano_data_significativos$Fosfopeptido, 
                               Tipo_Cancer = volcano_data_significativos$Tipo_Cancer)

# Encontrar modificaciones de secuencia comunes
fosfopeptidos_comunes <- intersect(mod_fosfopeptidos$Fosfopeptido, mod_volcano_plot$Fosfopeptido)

# Filtrar las filas de ambos data frames que contienen los fosfopeptidos comunes
mod_fosfopeptidos_comunes <- mod_fosfopeptidos[mod_fosfopeptidos$Fosfopeptido %in% fosfopeptidos_comunes, ]
mod_volcano_plot_comunes <- mod_volcano_plot[mod_volcano_plot$Fosfopeptido %in% fosfopeptidos_comunes, ]

# Unir las filas filtradas de ambos data frames
modificaciones_comunes_completas <- merge(mod_fosfopeptidos_comunes, mod_volcano_plot_comunes, 
                                          by = "Fosfopeptido", 
                                          suffixes = c("_fosfopeptidos", "_volcano"))

# Mostrar las modificaciones comunes combinadas
cat("\nModificaciones comunes combinadas:\n")
print(modificaciones_comunes_completas)

