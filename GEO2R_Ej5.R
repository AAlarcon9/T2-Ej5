#Ejercicio 5
#Geo2R
################################################################
#   Data plots for selected GEO samples
library(GEOquery)
library(limma)
library(umap)

# Cargar los datos de tu interés

gset <- getGEO("GSE68849", GSEMatrix =TRUE, getGPL=FALSE)

#Organizar los datos
if (length(gset) > 1) idx <- grep("GPL10558", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# Seleccipon de los nombres de las columnas 
fvarLabels(gset) <- make.names(fvarLabels(gset))

#Generar un vector de caracteres 
gsms <- "1010101010"

#Separar el vector de caracteres previamente definido
sml <- strsplit(gsms, split="")[[1]]

# Transformación log2 
#Recuperar los datos de expresión
ex <- exprs(gset)

qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))

LogC <- (qx[5] > 100) ||
  (qx[6]-qx[1] > 50 && qx[2] > 0)

if (LogC) { ex[which(ex <= 2)] <- NaN
exprs(gset) <- log2(ex) }

# Asignar las muestras a un grupo y crear la matriz
gs <- factor(sml)

#Grupos
groups <- make.names(c("Control","Infected"))

#Los niveles son los grupos
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)

#Ajuste a un modelo líneal
fit <- lmFit(gset, design)  # fit linear model
fit

# Indicar contrastes de interés and recalculate model coefficients
cts <- c(paste(groups[1],"-",groups[2],sep=""))
cont.matrix <- makeContrasts(contrasts=cts, levels=design)

# Volver a calcular los coeficientes del modelo
fit2 <- contrasts.fit(fit, cont.matrix)

# Obtener estadisticos
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=250)

#Recuperar tabla con los genes más signficativos

write.table(tT, file=stdout(), row.names=F, sep="\t")
write.table(tT,file="~/6to semestre/Genomica/Influenza.txt",sep="\t")

# Visualización (gráficas) ##########
# Elaborar histograma
tT2 <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
hist(tT2$adj.P.Val, col = "grey", border = "white", xlab = "P-adj",
     ylab = "Number of genes", main = "P-adj value distribution")

# separar los recultados como sobre, sub o no expresados
dT <- decideTests(fit2, adjust.method="fdr", p.value=0.01)

# Elaboración de un diagrama de Venn
vennDiagram(dT, circle.col=palette())

# Q-Q plot para obtener el t-statistic
#Filtrar los datos
t.good <- which(!is.na(fit2$F)) 
qqt(fit2$t[t.good], fit2$df.total[t.good], main="t-statistic")


# volcano plot
#Nombre de los grupos a comparar
colnames(fit2)

#indicar contraste
ct <- 1    

#Generar la gráfica
volcanoplot(fit2, coef=ct, main=colnames(fit2)[ct], pch=20,
            highlight=length(which(dT[,ct]!=0)), names=rep('+', nrow(fit2)))


################################################################