## ----setup, include=FALSE-----------------------------
#AnnArbor 
knitr::opts_chunk$set(echo = TRUE, comment = NA,cache=TRUE,out.width = "60%",fig.align='center')
leccion=0


## ---- echo=FALSE,fig.cap="",out.width = "20%"---------
knitr::include_graphics("Imgs/Rlogo.png")


## ----grafico_Rstudio2, echo=FALSE,fig.align='center',out.width = "20%"----
knitr::include_graphics("Imgs/RSLogo.png")


## ----rstudioID, echo=FALSE,out.width="55%",fig.align='center'----
knitr::include_graphics("Imgs/InterfazRStudio.png")


## ----icono_Rstudio3, echo=FALSE,fig.align='center',out.width="25%"----
knitr::include_graphics("Imgs/RSLogo.png")
knitr::include_graphics("Imgs/Rlogo.png")


## ----grafico_Rstudio_02, echo=FALSE,fig.align='center',out.width="60%"----
knitr::include_graphics("Imgs/easy_plus_tools.png")


## ----grafico_Rstudio, include=FALSE,fig.align='center'----
knitr::include_graphics("Imgs/help.png")


## -----------------------------------------------------
2+2
77%/%5
77%%5


## -----------------------------------------------------
sqrt(9)
log(exp(1))
log(1000,10)
log10(1000)


## -----------------------------------------------------
factorial(5)
choose(4,2)
factorial(6)
factorial(5)*6


## -----------------------------------------------------
sin(pi/2)
cos(pi)
tan(0)


## ----plot_tema1, eval=FALSE---------------------------
## x = seq(0,2*pi,0.1)
## plot(x,sin(x),type="l",col="blue",lwd=3,
##      xlab=expression(x), ylab="",
##      xlim=c(0,4),cex=0.5)
## curve(cos(x),col="red",add=TRUE)
## lines(x, tan(x/2), col="purple",lwd=3)
## legend("bottomleft",
##        col=c("blue","green","purple"),
##        legend=c("Seno","Coseno", "Tangente"),
##        lwd=3, bty="l",cex=0.8)


## ----plot2_tema1, fig.align = "center",echo=FALSE,out.width="60%"----
x = seq(0,2*pi,0.1)
plot(x,sin(x),type="l",col="blue",lwd=3, xlab=expression(x), ylab="",
     xlim=c(0,4),cex=0.5)
curve(cos(x),col="red",add=TRUE)
lines(x, tan(x/2), col="purple",lwd=3)
legend("bottomleft",col=c("blue","green","purple"),
     legend=c("Seno","Coseno", "Tangente"), lwd=3, bty="l",cex=0.8)


## -----------------------------------------------------
print(pi,5)
round(pi,5)
floor(pi)
ceiling(pi)


## -----------------------------------------------------
a= 8
cubo = function(x){x^3}
cubo(x=a)
raiz_cúbica = function(x){x^(1/3)}
raiz_cúbica(a)
raiz_cúbica(cubo(x=a))


## ----echo=FALSE,out.width="15%"-----------------------
knitr::include_graphics("Imgs/chunk_00.png")


## ---- comment = NA------------------------------------
x = 1+1
x


## ----echo=FALSE,out.width="30%",echo=FALSE------------
knitr::include_graphics("Imgs/primer_chunk.png")


## ----echo=FALSE,out.width="30%",echo=FALSE------------
knitr::include_graphics("Imgs/segundo_chunk.png")


## ---- echo=FALSE,out.width="30%"----------------------
knitr::include_graphics("Imgs/no_aparece.png")


## ---- echo =FALSE-------------------------------------
sec = 10:20
sec
cumsum(sec)


## ---- echo=FALSE, message = TRUE,out.width="30%",echo=FALSE----
knitr::include_graphics("Imgs/parametros_chunk_2.png")


## ---- echo = TRUE, message = TRUE---------------------
library(car)
head(cars,3)


## ---- echo=FALSE, out.width="30%",echo=FALSE----------
knitr::include_graphics("Imgs/para_chunks_3.png")


## ---- echo = TRUE, message = FALSE, comment = NA------
library(car)
head(cars,3)


## ----out.width="30%",echo=FALSE-----------------------
knitr::include_graphics("Imgs/para_chunks_03.png")


## ---- echo=TRUE, results='markup'---------------------
sec = 10:20
sec
cumsum(sec)


## ----out.width="30%",echo=FALSE-----------------------
knitr::include_graphics("Imgs/parametros_chunk_4.png")


## ---- echo=TRUE, results='hide'-----------------------
sec = 10:20
sec
cumsum(sec)


## ----out.width="30%",echo=FALSE-----------------------
knitr::include_graphics("Imgs/parametros_chunk_5.png")


## ----chunk_ex, echo=TRUE, results='asis'--------------
sec = 10:20
sec
cumsum(sec)


## ----out.width="30%",echo=FALSE-----------------------
knitr::include_graphics("Imgs/parametros_chunk_6.png")


## ----una_chunk, echo=TRUE, results='hold'-------------
sec = 10:20
sec
cumsum(sec)


## -----------------------------------------------------
c(1,2,3)
rep("Mates",7)


## ----scan3_copias,out.width="60%",echo=FALSE----------
knitr::include_graphics("Imgs/scan.png")


## -----------------------------------------------------
cuadrado = function(x){x^2}
v = c(1,2,3,4,5,6)
sapply(v, FUN = cuadrado)
mean(v)
cumsum(v)


## -----------------------------------------------------
v = c(1,7,5,2,4,6,3)
sort(v)
rev(v)


## -----------------------------------------------------
v = c(14,5,6,19,32,0,8)
v[2]
v[-c(3,5)]
v[v != 19 & v>15]


## ----echo=FALSE---------------------------------------
cuadrado = function(x){x^2}
v = c(1,2,3,4,5,6)
sapply(v, FUN = cuadrado)
mean(v)
cumsum(v)


## -----------------------------------------------------
fac = factor(c(1,1,1,2,2,3,2,4,1,3,3,4,2,3,4,4), 
       levels = c(1,2,3,4), 
       labels = c("Sus","Apr","Not","Exc"))
fac
facOrd = ordered(c(1,1,1,2,2,3,2,4,1,3,3,4,2,3,4,4), 
       levels = c(1,2,3,4), 
       labels = c("Sus","Apr","Not","Exc"))
facOrd


## -----------------------------------------------------
x = c(1,-2,3,4,-5,6,7,-8,-9,0)
miLista = list(nombre = "X", vector = x, media = mean(x), 
               sumas = cumsum(x))
miLista


## -----------------------------------------------------
str(miLista)
names(miLista)


## -----------------------------------------------------
A = matrix(c(1,2,3,4,5,6,7,8,9), ncol = 3)
dim(A)
diag(A)


## -----------------------------------------------------
apply(A, MARGIN = c(1,2), FUN = cuadrado)
apply(A, MARGIN = 1, FUN = sum)
apply(A, MARGIN = 2, FUN = sum)


## -----------------------------------------------------
M = rbind(c(2,6,-8), c(0,6,-3), c(0,2,1))
eigen(M)


## -----------------------------------------------------
M = matrix(c(0,1,0,-7,3,-1,16,-3,4), nrow=3, byrow=TRUE)
eigen(M)


## ---- fig.height = 4, fig.width = 7, fig.align = "center",out.width="50%"----
alumnos = c(1:10)
notas = c(2,5,7,9,8,3,5,6,10,7)
plot(alumnos,notas)


## ----cph , echo=FALSE,fig.align = "center", out.width="45%"----
knitr::include_graphics("Imgs/pch.png")


## ---- fig.height = 3.75, fig.width = 9, fig.align = "center",out.width="60%"----
par(mfrow = c(1,2))
plot = plot(exp(1:20), xlab = "Indice",
            ylab = expression(e^{1:20}), 
            main = "Escala lineal")
plotLog = plot(exp(1:20), log = "y", xlab = "Indice", 
               ylab = expression(e^{1:20}), 
               main = "Escala logaritmica en el eje y")
par(mfrow = c(1,1))


## ---- eval=FALSE--------------------------------------
## par(mfrow = c(3,2))
## x = c(50:59)
## y = c(2,9,25,3,100,77,62,54,19,40)
## plot(x,y, pch = 23, cex = 2, col = "blue", type = "p")
## plot(x,y, pch = 23, cex = 2, col = "blueviolet", type = "l")
## plot(x,y, pch = 23, cex = 2, col = "gold", type = "b")
## plot(x,y, pch = 23, cex = 2, col = "deeppink", type = "o")
## plot(x,y, pch = 23, cex = 2, col = "springgreen",
##      type = "h")
## plot(x,y, pch = 23, cex = 2, col = "firebrick1",
##      type = "s")
## par(mfrow = c(1,1))


## ----tipo2,fig.align='center',out.width="45%",echo=FALSE----
knitr::include_graphics("Imgs/tipo_grafico_1.png")


## ---- fig.height = 4, fig.width = 9, fig.align = "center",out.width="60%"----
x = (2*(1:20))
y = (-1)^(1:20)*5*(1:20)
plot(x,y, main = "Ejemplo de grafico", pch = 8, cex = 1,
     type = "b", lty = 4, lwd = 4, 
     xaxp = c(0,40,2), yaxp = c(-100,100,8))


## ---- eval=FALSE--------------------------------------
## x = (2*(1:20))
## y = (-1)^(1:20)*5*(1:20)
## plot(x,y, main = "Poniendo un punto y una recta", pch = 8,
##      cex = 1, type = "b", lty = 4,
##      lwd = 4, xaxp = c(0,40,2), yaxp = c(-100,100,8))
## points(20,0, col = "red", cex = 4, pch = 16)
## abline (h = 0, lty = 2, col = "dodgerblue")


## ---- echo = FALSE, fig.width = 9, fig.align = "center",out.width="50%"----
x = (2*(1:20))
y = (-1)^(1:20)*5*(1:20)
plot(x,y, main = "Poniendo un punto y una recta", pch = 8, 
     cex = 1, type = "b", lty = 4, lwd = 4, 
     xaxp = c(0,40,2), yaxp = c(-100,100,8))
points(20,0, col = "red", cex = 4, pch = 16)
abline (h = 0, lty = 2, col = "dodgerblue")


## ---- fig.width = 9, fig.height = 3.75,fig.align = "center",out.width="65%"----
alumnos = c(1:10)
notas = c(2,5,7,9,8,3,5,6,10,7)
plot(alumnos,notas, main = "Grafico con texto")
text(alumnos,notas, 
     labels = c("S","A","N","E","N","S","A","A","E","N"), 
     pos = c(rep(3,times = 8),1,3))


## ---- results='hide', fig.align="center", fig.height=4,out.width="80%"----
x = c(5*(1:20))
plot(x,c(exp(-x)+(-1)^x*x/2*sin(x)^2))
lines(c(20,10,40,80,60,60,20),c(20,0,-20,-20,40,0,20), 
      lwd = 2, col = "darkslategray1")
curve(20*sin(x), add = TRUE, col = "green")


## ---- eval = FALSE------------------------------------
## x = seq(0,2*pi,0.1)
## plot(x,sin(x),type="l",col="blue",lwd=3, xlab="", ylab="")
## lines(x,cos(x),col="green",lwd=3)
## lines(x, tan(x), col="purple",lwd=3)
## legend("bottomleft",col=c("blue","green","purple"),
##        legend=c("Seno","Coseno", "Tangente"),
##        lwd=3, bty="l")


## ---- echo = FALSE, fig.align="center",out.width="60%"----
x = seq(0,2*pi,0.1)
plot(x,sin(x),type="l",col="blue",lwd=3,
     xlab="", ylab="")
lines(x,cos(x),col="green",lwd=3)
lines(x, tan(x), col="purple",lwd=3)
legend("bottomleft",
       col=c("blue","green","purple"),
       legend=c("Seno","Coseno", "Tangente"),
       lwd=3, bty="l")


## ---- eval = FALSE------------------------------------
## x = c(5*(1:10))
## plot(x,c(exp(-x)+(-1)^x*x/2*sin(x)^2), xlab = "", ylab = "",
##      main = "Grafico con varios elementos")
## segments(10,0,40,0, col = "red", lwd = 4)
## arrows(10,0,40,-10, col = " blue", length = 0.5,
##        angle = 5, code = 3)
## symbols(40,0,stars = cbind(1,.5,1,.5,1,.5,1,.5,1,.5),
##         add = TRUE, lwd = 3, inches = 0.5)
## symbols(40,0,stars = cbind(1,.5,1,.5,1,.5,1,.5,1,.5),
##         add = TRUE, lwd = 3)
## polygon(c(20,30,40),c(10,-10,10), col = "gold",
##         density = 3, angle = 90, lty = 4,
##         lwd = 5)


## ---- echo = FALSE, fig.align="center",out.width="60%"----
x = c(5*(1:10))
plot(x,c(exp(-x)+(-1)^x*x/2*sin(x)^2), xlab = "", ylab = "", 
     main = "Grafico con varios elementos")
segments(10,0,40,0, col = "red", lwd = 4)
arrows(10,0,40,-10, col = " blue", length = 0.5, angle = 5, code = 3)
symbols(40,0,stars = cbind(1,.5,1,.5,1,.5,1,.5,1,.5),
        add = TRUE, lwd = 3, inches = 0.5)
symbols(40,0,stars = cbind(1,.5,1,.5,1,.5,1,.5,1,.5),
        add = TRUE, lwd = 3)
polygon(c(20,30,40),c(10,-10,10), col = "gold", density = 3, 
        angle = 90, lty = 4, lwd = 5)


## -----------------------------------------------------
str(Orange)


## -----------------------------------------------------
head(Orange,4)
tail(Orange,4)


## -----------------------------------------------------
dataOrange = Orange
dataOrange[c(10:12),]
dataOrange[c(2,17),c(1,3)]


## -----------------------------------------------------
dataOrange[2,3]
dataOrange[dataOrange$circumference<=50,]


## -----------------------------------------------------
notas = read.table(
  "http://aprender.uib.es/Rdir/Controls11-12.txt", 
  col.names = c("Nota_Parcial","Nota_Final","Grup"),
  sep=",",header=TRUE)
head(notas,8)


## -----------------------------------------------------
write.table(notas, file = "data/NotasData.csv",
            dec = ".")
notas2 = read.table("data/NotasData.csv", header = TRUE)
str(notas2)


## -----------------------------------------------------
Programacion = c(1,2,0,5,4,6,7,5,5,8)
Calculo = c(3,3,2,7,9,5,6,8,5,6)
Empresa = c(4,5,4,8,8,9,6,7,9,10)
grados = data.frame(Pr = Programacion, 
                    Ca = Calculo, Em = Empresa)
str(grados)


## -----------------------------------------------------
Ingles = c(5,4,6,2,1,0,7,8,9,6)
grados2 = cbind(grados, Ingles)
head(grados2)


## -----------------------------------------------------
head(iris ,5)
str(iris)


## -----------------------------------------------------
x = sample(1:5, size = 12, replace = TRUE)
x

Respuestas=factor(sample(c("Si", "No"), size = 12, replace = TRUE)) 
Respuestas


## -----------------------------------------------------
table(x)

table(Respuestas)


## -----------------------------------------------------
names(table(x))

names(table(Respuestas))


## -----------------------------------------------------
z=factor(x, levels=1:7) #Los niveles serán 1,2,3,4,5,6,7 
z
table(z)


## -----------------------------------------------------
table(x)[3] #La tercera columna de table(x)
table(x)["7"] #¿La columna de table(x) con nombre 7?


## -----------------------------------------------------
table(x)["5"] #La columna de table(x) con nombre 5
3*table(x)[2] #El triple de la segunda columna de table(x)


## -----------------------------------------------------
sum(table(x)) #Suma de las entradas de table(x)
sqrt(table(Respuestas))
# Raíces cuadradas de las entradas de table(Respuestas)


## -----------------------------------------------------
prop.table(table(x))

prop.table(table(Respuestas))


## -----------------------------------------------------
prop.table(x)


## -----------------------------------------------------
X=c(1,1,1)
prop.table(table(X))
prop.table(X)


## -----------------------------------------------------
table(x)/length(x)



## -----------------------------------------------------
table(x)
names(which(table(x)==1))


## -----------------------------------------------------
names(which(table(x)==max(table(x))))
names(which(table(Respuestas)==max(table(Respuestas))))


## -----------------------------------------------------
Sexo_Ger=c("Mujer","Mujer","Hombre","Mujer","Mujer","Mujer",
           "Mujer","Mujer","Hombre","Mujer","Hombre","Hombre",
           "Mujer", "Mujer","Hombre","Mujer","Mujer","Mujer",
           "Mujer","Hombre")
t0=table(Sexo_Ger)
prop.table(t0) 
names(which(t0==max(t0))) 


## -----------------------------------------------------
Sexo= sample(c("H", "M"), size = length(Respuestas), 
             replace = T) #H = hombre, M = mujer
table(Respuestas ,Sexo)


## -----------------------------------------------------
table(Respuestas ,Sexo)[1,2]
table(Respuestas ,Sexo)["No","M"]


## -----------------------------------------------------
prop.table(table(Sexo,Respuestas)) #Global


## -----------------------------------------------------
prop.table(table(Sexo,Respuestas), margin=1) #Por sexo
prop.table(table(Sexo,Respuestas), margin=2) #Por respuesta


## -----------------------------------------------------
table(Sexo,Respuestas) 


## -----------------------------------------------------
colSums(table(Sexo,Respuestas)) 
rowSums(table(Sexo,Respuestas)) 


## -----------------------------------------------------
colSums(prop.table(table(Sexo,Respuestas)))
rowSums(prop.table(table(Sexo,Respuestas)))


## -----------------------------------------------------
Beb_Energ=read.table("data/EnergyDrink",header=TRUE)


## -----------------------------------------------------
str(Beb_Energ)
head(Beb_Energ,4)


## -----------------------------------------------------
summary(Beb_Energ)


## -----------------------------------------------------
summary(Beb_Energ)[,2]


## -----------------------------------------------------
apply(Beb_Energ, MARGIN=2, FUN=table)


## -----------------------------------------------------
apply(Beb_Energ,MARGIN=2,FUN=table)$sexo

table(Beb_Energ$sexo)


## -----------------------------------------------------
table(Beb_Energ)


## -----------------------------------------------------
table(Beb_Energ[c(1,3)])


## -----------------------------------------------------
ftable(Beb_Energ)


## ----out.width="50%",fig.align='center',eval=FALSE----
## barplot(table(Sexo), col=c("lightblue","pink"),
## main="Diagrama de barras de
## las frecuencias absolutas\n de la variable \"Sexo\"")


## ----out.width="45%",fig.align='center',echo=FALSE----
barplot(table(Sexo), col=c("lightblue","pink"), 
main="Diagrama de barras de 
las frecuencias absolutas\n de la variable \"Sexo\"")


## ----out.width="50%",fig.align='center',eval=FALSE----
## barplot(prop.table(table(Respuestas)),
## main="Diagrama de barras de frecuencias
##       relativas\n de la variable \"Respuestas\"")


## ----out.width="45%",fig.align='center',echo=FALSE----
barplot(prop.table(table(Respuestas)), 
main="Diagrama de barras de frecuencias 
      relativas\n de la variable \"Respuestas\"")


## ----out.width="45%",fig.align='center'---------------
par(mfrow=c(1,2))
barplot(table(Respuestas), col=c("green"))
barplot(table(Respuestas), col=c("red","blue"))
par(mfrow=c(1,1))


## ----out.width="50%",fig.align='center'---------------
barplot(table(x), horiz=TRUE)


## ----out.width="50%",fig.align='center'---------------
barplot(table(Sexo,Respuestas), legend.text = TRUE)


## ----out.width="50%",fig.align='center'---------------
barplot(table(Sexo,Respuestas), beside=TRUE, 
        legend.text=TRUE)


## ----out.width="50%",fig.align='center'---------------
barplot(table(Respuestas,Sexo), beside=TRUE, 
        names=c("Men", "Women"), col=c("yellow","lightblue"),  legend.text=c("No","Yes"))


## ----out.width="20%",fig.align='center'---------------
pie(table(x), main="Diagrama circular de la variable x")


## ----out.width="20%",fig.align='center'---------------
Respuestas
pie(table(Respuestas), main="Diagrama circular de la variable Respuestas")


## ----out.width="40%",fig.align='center'---------------
plot(table(Sexo,Respuestas), 
main="Gráfico de mosaico de las variables
  \"Sexo\" y \"Respuestas\"")


## ----out.width="40%",fig.align='center'---------------
plot(HairEyeColor, main="Gráfico de mosaico de la tabla HairEyeColor", 
     col=c("pink","lightblue"))


## -----------------------------------------------------
#library(vcd)
#País = sample(c("Francia","Alemania","España"), size = length(Sexo), replace = T)
#cotabplot(table(Sexo,Respuestas,País))


## -----------------------------------------------------
#library(vcdExtra)
#mosaic3d(HairEyeColor, type="expected", box=TRUE,
#col=c("pink","lightblue"))


## ---- echo=FALSE--------------------------------------
HEC=as.table(HairEyeColor[ , , 1]+ HairEyeColor[ , , 2])


## -----------------------------------------------------
dimnames(HEC)


## ---- echo = FALSE------------------------------------
dimnames(HEC)=list(Cabello=c("Negro","Marron","Rojo","Rubio"), Ojos=c("Marrones","Azules","Pardos","Verdes"))


## ----mosaico_x1 ,out.width="50%",fig.align='center',echo = FALSE----
plot(HEC,col=c("lightblue"),main="Diagrama de mosaico de la tabla bidimensional de frecuencias
     \n de colores de cabellos y ojos")


## ---- echo = FALSE------------------------------------
sum(HEC) #Número total de individuos [1] 592


## ---- echo = FALSE------------------------------------
colSums(HEC) #Frec. abs. de Ojos
rowSums(HEC)
round(prop.table(colSums(HEC)),3) #Frec. rel. de Ojos
round(prop.table(rowSums(HEC)),3) #Frec. rel. de Cabello Negro Castaño Rojo Rubio


## ----out.width="50%",fig.align='center', echo = F-----
par(mfrow=c(1,2))
barplot(prop.table(colSums(HEC)), ylim=c(0,0.4), col=c("burlywood4","lightblue","orange3","lightgreen"), main="Frecuencias relativas \nde colores de ojos")
barplot(prop.table(rowSums(HEC)),
        col=c("black","chocolate4","red","gold"), ylim=c(0,0.5), main="Frecuencias relativas \nde colores de cabello")
par(mfrow=c(1,1))


## ---- echo = F----------------------------------------
round(prop.table(HEC), 3) #Frec. rel. globales


## ---- echo = F----------------------------------------

round(prop.table(HEC, margin=1), 3) #Frec. rel. por filas
round(prop.table(HEC, margin=2), 3) #Frec. rel. por columnas


## ---- echo = F,out.width="50%",fig.align='center'-----
par(mfrow=c(1,2))
barplot(prop.table(HEC, margin=1), beside=TRUE, legend.text=TRUE,
col=c("black","brown","red","gold"), ylim=c(0,0.8), main="Frecuencias relativas de colores de \ncabello en cada color de ojos")

barplot(t(prop.table(HEC, margin=2)), beside=TRUE, legend.text=TRUE, ylim=c(0,0.6), col=c("burlywood4","lightblue","orange3","lightgreen"), main="Frecuencias relativas de colores \nde ojo en cada color de cabellos")
par(mfrow = c(1,1))


## -----------------------------------------------------
notas = ordered(c("S","A", "N", "Ex", "S", "S",
                  "Ex","Ex", "N", "A", "A", "A",
                  "A", "N", "S"),
                levels = c("S", "A", "N", "Ex"))
table(notas)


## -----------------------------------------------------
set.seed(2018)
clientes = sample(1:5, 50, replace = TRUE)
clientes
set.seed(NULL)


## ---- echo = FALSE------------------------------------
absolut = table(clientes)
relative = prop.table(absolut)
acumul = cumsum(absolut)
rel.acumul = cumsum(relative)
absolut = (as.matrix(absolut))
relative = (as.matrix(relative))
acumul = (as.matrix(acumul))
rel.acumul = (as.matrix(rel.acumul))

clientela = data.frame(absolut,relative,acumul,rel.acumul)
colnames(clientela) = c("Absoluta", "Relativa", "Acumulada", "Rel. Acumulada")
clientela



## -----------------------------------------------------
notas
fAbs = table(notas) #Frec. abs.
cumsum(fAbs) #Frec. abs. acumuladas


## ---- fig.height=3.75,,out.width="50%",fig.align='center'----
cumsum(prop.table(fAbs)) #Frec. relativas acumuladas
barplot(fAbs, main = "Diagrama de barras de frecuencias absolutas")


## ----out.width="50%",fig.align='center'---------------
barplot(cumsum(fAbs), main = "Diagrama de barras de 
    frecuencias absolutas acumuladas")


## -----------------------------------------------------
cumsum(table(notas))/length(notas)
cumsum(table(notas)/length(notas))


## ---- echo = FALSE------------------------------------
set.seed(2018)
longitud = sample(1:5,100, replace = TRUE)
longitud = ordered(longitud)
levels(longitud) = c("Muy.corto","Corto","Normal","Largo","Muy.largo")


## -----------------------------------------------------
longitud

## ---- echo = FALSE------------------------------------
set.seed(NULL)


## -----------------------------------------------------
Fr.Abs = table(longitud)
Fr.Abs
Fr.Rel = prop.table(Fr.Abs)
Fr.Rel


## -----------------------------------------------------
Fr.Acum = cumsum(Fr.Abs)
Fr.Acum
Fr.RAcum = cumsum(Fr.Rel)
Fr.RAcum


## ----out.width="40%",fig.align='center'---------------
barplot(Fr.RAcum, main = "Diagrama de frecuencias relativas acumuladas")


## -----------------------------------------------------
zonas = rep(c("A","B","C","D"), c(30,25,35,10))
jirafas = data.frame(zonas,longitud)
str(jirafas)
head(jirafas)


## -----------------------------------------------------
apply(table(jirafas), MARGIN = 1, FUN = cumsum)


## -----------------------------------------------------
t(apply(table(jirafas), MARGIN = 1, FUN = cumsum))


## -----------------------------------------------------
t(apply(prop.table(table(jirafas),
                   margin = 1), 
        MARGIN = 1, FUN = cumsum))


## ----eval=FALSE,,out.width="50%",fig.align='center'----
## Diagrama = apply(prop.table(table(jirafas), margin = 1),
##                  MARGIN = 1, FUN = cumsum)
## barplot(Diagrama, beside = TRUE, legend = TRUE,
## main = "Diagrama de barras de frecuencias relativas acumuladas
##         de longitudes por zonas",
## args.legend=list(x="topleft", cex=0.55))


## ---- echo=FALSE,out.width="50%",fig.align='center'----
Diagrama = apply(prop.table(table(jirafas), margin = 1),
                 MARGIN = 1, FUN = cumsum)
barplot(Diagrama, beside = TRUE, legend = TRUE,
main = "Diagrama de barras de frecuencias relativas acumuladas 
        de longitudes por zonas",
args.legend=list(x="topleft", cex=0.55))


## -----------------------------------------------------
crabs = read.table("data/datacrab.txt", header = TRUE)
crabs = crabs[,-1] #Omitimos la primera columna
str(crabs)


## -----------------------------------------------------
table(crabs$width)


## -----------------------------------------------------
intervalos = cut(crabs$width, breaks = c(21,25,29,33,Inf), right = FALSE, 
                 labels = c("21-25", "25-29", "29-33", "33-..."))


## -----------------------------------------------------
crabs$width.rank = ordered(intervalos)
str(crabs)


## -----------------------------------------------------
Tabla = table(crabs[,c(1,6)])
Tabla


## -----------------------------------------------------
Fr.rel = round(prop.table(Tabla,margin = 1),3)
Fr.rel


## -----------------------------------------------------
Fr.rel.acu = round(apply(prop.table(Tabla, margin = 1), 
                         MARGIN = 1, FUN = cumsum), 3)
t(Fr.rel.acu)


## ---- eval=FALSE--------------------------------------
## azul = c("cyan", "cyan4", "cyan1", "cyan3")
## 
## barplot(t(Fr.rel), beside = TRUE, legend = TRUE, ylim = c(0,1), col = azul,
##         main = "Diagrama de barras de frecuencias relativas",
##         args.legend=list(x = "topright", cex=0.55))


## ----echo=FALSE ,out.width="50%",fig.align='center'----
azul = c("cyan", "cyan4", "cyan1", "cyan3")

barplot(t(Fr.rel), beside = TRUE, legend = TRUE, ylim = c(0,1), col = azul, 
        main = "Diagrama de barras de frecuencias relativas", 
        args.legend=list(x = "topright", cex=0.55))


## ---- eval=FALSE--------------------------------------
## barplot(Fr.rel.acu, beside = TRUE, legend = TRUE, col = azul,
##         main = "Diagrama de barras de frecuencias relativas acumuladas",
##         args.legend=list(x = "topleft", cex=0.55))


## ----out.width="50%",fig.align='center',echo=FALSE----
barplot(Fr.rel.acu, beside = TRUE, legend = TRUE, col = azul, 
        main = "Diagrama de barras de frecuencias relativas acumuladas", 
        args.legend=list(x = "topleft", cex=0.55))


## -----------------------------------------------------
edad = c(15,18,25,40,30,29,56,40,13,27,42,23,11,26,25,32,30,40,33,29)


## -----------------------------------------------------
table(edad)


## -----------------------------------------------------
round(prop.table(table(edad)),3)
cumsum(table(edad))


## -----------------------------------------------------
round(cumsum(prop.table(table(edad))),3)


## -----------------------------------------------------
set.seed(162017)
dados = sample(1:6,25,replace = TRUE)
dados
set.seed(NULL)


## -----------------------------------------------------
table(dados)
round(prop.table(table(dados)),2)
cumsum(table(dados))


## -----------------------------------------------------
round(cumsum(prop.table(table(dados))),2)
dados.df = data.frame(
  Puntuacion = 1:6,
  Fr.abs = as.vector(table(dados)),
  Fr.rel = as.vector(round(prop.table(table(dados)),2)),
  Fr.acu = as.vector(cumsum(table(dados))),
  Fr.racu = as.vector(round(cumsum(
    prop.table(table(dados))),2)))


## -----------------------------------------------------
dados.df


## -----------------------------------------------------
sort(edad) #Ordenamos los datos por su orden natural
table(edad)


## -----------------------------------------------------
dados.df


## -----------------------------------------------------
mean(edad) #La media aritmética
mean(dados)
median(edad) #La mediana


## -----------------------------------------------------
median(dados)
as.numeric(names(which(
  table(edad)==max(table(edad))))) #La moda
as.numeric(names(which(
  table(dados)==max(table(dados)))))


## -----------------------------------------------------
set.seed(260798)
dado = sample(1:4, 50, replace = TRUE)
set.seed(NULL)
length(dado)
dado = sort(dado) #Los ordenamos de menor a mayor
dado


## -----------------------------------------------------
df.dado = data.frame(
  Puntuacion = 1:4,
  Fr.abs = as.vector(table(dado)),
  Fr.rel = as.vector(round(prop.table(table(dado)),2)),
  Fr.acu = as.vector(cumsum(table(dado))),
  Fr.racu = as.vector(round(cumsum(
    prop.table(table(dado))),2))
  )
df.dado


## -----------------------------------------------------
dado[15]


## -----------------------------------------------------
set.seed(0)
dados2 = sample(1:6,15, replace = TRUE)
dados2
set.seed(NULL)
quantile(dados2,0.25) #Primer cuartil
quantile(dados2,0.8)


## ----echo=FALSE---------------------------------------
set.seed(260798)
dado = sample(1:4, 50, replace = TRUE)
set.seed(NULL)
dado = sort(dado) #Los ordenamos de menor a mayor
set.seed(162017)
dados = sample(1:6,25,replace = TRUE)
set.seed(NULL)
set.seed(0)
dados2 = sample(1:6,15, replace = TRUE)
set.seed(NULL)


## -----------------------------------------------------
dados2
diff(range(dados2))
IQR(dados2)
var(dados2)


## -----------------------------------------------------
sd(dados2)
n = length(dados2)
var(dados2)*(n-1)/n
sd(dados2)*sqrt((n-1)/n)


## -----------------------------------------------------
cangrejos = read.table("data/datacrab.txt", header = TRUE) 
#Cargamos el data frame
cangrejos = cangrejos[-1] #Eliminamos la primera columna
summary(cangrejos) #Aplicamos la función summary


## -----------------------------------------------------
summary(subset(cangrejos, color == 3,c("weight","width")))



## -----------------------------------------------------
summary(subset(cangrejos, color == 5,c("weight","width")))


## ---- results="hide"----------------------------------
by(iris[,c(1,3)], iris$Species, FUN = summary)


## ---- results="hide"----------------------------------
aggregate(cbind(Sepal.Length,Petal.Length)~Species, data=iris, FUN=summary)


## -----------------------------------------------------
dadosNA = c(dados2,NA)
dadosNA
mean(dadosNA)
mean(dadosNA, na.rm = TRUE)


## ----out.width="45%",fig.align='center'---------------
boxplot(dados2, main = "Un diagrama de caja")


## ----out.width="50%",fig.align='center'---------------
boxplot(dado,dados,dados2)


## ---- fig.width=10, fig.height=5, out.width="80%",fig.align='center'----
body = read.table("data/bodyfat.txt", header = TRUE)
boxplot(body,las=2)


## ----out.width="50%",fig.align='center'---------------
boxplot(body[,7:9], names = c("Pecho", "Abdomen", "Cadera"))


## ----out.width="50%",fig.align='center'---------------
boxplot(circumference~Tree, data = Orange, ylab = "Circunferencia del tronco (mm)", 
        main = "Boxplot de los naranjos en función del tipo de árbol")


## ----out.width="50%",fig.align='center'---------------
boxplot(Sepal.Width~Species, data = iris, ylab = "Anchura del sétalo (cm)",
        notch = TRUE, col = c("cyan","cyan2","cyan4"),
        main = "Boxplot de iris")


## ---- fig.height=3.5 , out.width="80%"----------------
boxplot(Sepal.Width~Species, data = iris, ylab = "Anchura del sétalo (cm)")
medias = aggregate(Sepal.Width~Species, data = iris, FUN = mean)
points(medias, col = "pink", pch = 15)


## ----echo=FALSE---------------------------------------
set.seed(260798)
dado = sample(1:4, 50, replace = TRUE)
set.seed(NULL)
dado = sort(dado) #Los ordenamos de menor a mayor
set.seed(162017)
dados = sample(1:6,25,replace = TRUE)
set.seed(NULL)
set.seed(0)
dados2 = sample(1:6,15, replace = TRUE)
set.seed(NULL)


## -----------------------------------------------------
pesos = c(55.2,54.0,55.2,53.7,60.2,53.2,54.6,55.1,51.2,53.2,54.8,52.3,
          56.9,57.0,55.0,53.5,50.9,55.1,53.6,61.2,59.5,50.3,52.7,60.0)


## ----out.width="50%",fig.align='center'---------------
barplot(table(pesos))


## ----eval=FALSE, out.width="60%",fig.align='center'----
## hist(pesos, breaks = 6, right = FALSE, main = "", ylab = "Frecuencia", xlab = "pesos")


## ----echo=FALSE ,out.width="50%",fig.align='center'----
hist(pesos, breaks = 6, right = FALSE, main = "", ylab = "Frecuencia", 
     xlab = "pesos")


## -----------------------------------------------------
crabs = read.table("data/datacrab.txt", header = TRUE)
str(crabs)
cw = crabs$width


## -----------------------------------------------------
n = length(cw)
k1 = ceiling(sqrt(n))
k1


## -----------------------------------------------------
k2 = ceiling(1+log(n,2))
k2


## -----------------------------------------------------
As = 3.5*sd(cw)*n^(-1/3) #Amplitud teórica
k3 = ceiling(diff(range(cw))/As)
k3


## -----------------------------------------------------
#Amplitud teórica
Afd = 2*(quantile(cw,0.75, names = FALSE)-
           quantile(cw,0.25,names = FALSE))*n^(-1/3) 
k4 = ceiling(diff(range(cw))/Afd)
k4


## -----------------------------------------------------
nclass.Sturges(cw)
nclass.scott(cw)
nclass.FD(cw)


## -----------------------------------------------------
A = diff(range(cw)) / 10
A


## -----------------------------------------------------
A = 1.3


## -----------------------------------------------------
L1 = min(cw)-1/2*0.1
L1


## -----------------------------------------------------
L2 = L1 + A
L3 = L2 + A
L4 = L3 + A
L5 = L4 + A
L6 = L5 + A
L7 = L6 + A
L8 = L7 + A
L9 = L8 + A
L10 = L9 + A
L11 = L10 + A
L = c(L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11)
L


## -----------------------------------------------------
L = L1 + A*(0:10)
L


## -----------------------------------------------------
X1 = (L[1]+L[2])/2
X1


## -----------------------------------------------------
X2 = X1 + A
X3 = X2 + A
X4 = X3 + A
X5 = X4 + A
X6 = X5 + A
X7 = X6 + A
X8 = X7 + A
X9 = X8 + A
X10 = X9 + A
X = c(X1,X2,X3,X4,X5,X6,X7,X8,X9,X10)
X


## -----------------------------------------------------
X = X1 + A*(0:9)
X


## -----------------------------------------------------
X = (L[1:length(L)-1]+L[2:length(L)])/2
X


## -----------------------------------------------------
#Primera función
TablaFrecs = function(x,k,A,p){ 
  L = min(x)-p/2+A*(0:k)
  x_cut = cut(x, breaks = L, right=FALSE)
  intervals = levels(x_cut)
  mc = (L[1]+L[2])/2+A*(0:(k-1))
  Fr.abs = as.vector(table(x_cut)) 
  Fr.rel = round(Fr.abs/length(x),4) 
  Fr.cum.abs = cumsum(Fr.abs) 
  Fr.cum.rel = cumsum(Fr.rel)
  tabla = data.frame(intervals, mc, Fr.abs, Fr.cum.abs, Fr.rel, Fr.cum.rel)
  tabla
  }


## -----------------------------------------------------
TablaFrecs.L = function(x,L,V){
  x_cut = cut(x, breaks=L, right=FALSE, include.lowest=V)
  intervals = levels(x_cut)
  mc = (L[1:(length(L)-1)]+L[2:length(L)])/2
  Fr.abs = as.vector(table(x_cut)) 
  Fr.rel = round(Fr.abs/length(x),4)
  Fr.cum.abs = cumsum(Fr.abs)
  Fr.cum.rel = cumsum(Fr.rel)
  tabla = data.frame(intervals, mc, Fr.abs, Fr.cum.abs, Fr.rel, Fr.cum.rel)
  tabla
  }


## ---- echo = FALSE------------------------------------
n = c(length(which(cw>=L[1] & cw<L[2])),length(which(cw>=L[2] & cw<L[3])),length(which(cw>=L[3] & cw<L[4])),
      length(which(cw>=L[4] & cw<L[5])),
      length(which(cw>=L[5] & cw<L[6])),
      length(which(cw>=L[6] & cw<L[7])),
      length(which(cw>=L[7] & cw<L[8])),
      length(which(cw>=L[8] & cw<L[9])),
      length(which(cw>=L[9] & cw<L[10])),
      length(which(cw>=L[10] & cw<L[11])))
N = cumsum(n)
f = round(n/N[10],4)
FF = cumsum(f)


## -----------------------------------------------------
TablaFrecs(cw,10,1.3,0.1)


## -----------------------------------------------------
TablaFrecs.L(cw,L,FALSE)


## ---- echo = FALSE------------------------------------
set.seed(4)
notas = sample(0:10,100, replace = TRUE)
set.seed(NULL)


## -----------------------------------------------------
notas


## -----------------------------------------------------
#Definimos vector de extremos
L = c(0,5,7,9,10)
#Definimos notas1 como el resultado de la codificación en 
#intervalos utilizando como etiquetas los propios intervalos
notas1 = cut(notas, breaks = L, right = FALSE, include.lowest = TRUE)
notas1


## -----------------------------------------------------
#Definimos las marcas de clase
MC = (L[1:length(L)-1]+L[2:length(L)])/2
#Definimos notas2 como el resultado de la codificación en 
#intervalos utilizando como etiquetas las marcas de clase
notas2 = cut(notas, breaks = L, labels = MC, right = FALSE, 
             include.lowest = TRUE)
notas2


## -----------------------------------------------------
#Definimos notas3 como el resultado de la codificación en 
#intervalos utilizando como etiquetas la posición ordenada 
#del intervalo (1, 2, 3 o 4)
notas3 = cut(notas, breaks = L, labels = FALSE, right = FALSE, 
             include.lowest = TRUE)
notas3


## -----------------------------------------------------
#Definimos notas4 como el resultado de la codificación en 
#intervalos utilizando como etiquetas Susp, Aprob, Not y Exc
notas4 = cut(notas, breaks = L, labels = c("Susp", "Aprob", "Not", "Exc"), 
             right = FALSE, include.lowest = TRUE)
notas4


## -----------------------------------------------------
cut(notas, breaks = 4, right = FALSE, include.lowest = TRUE)


## -----------------------------------------------------
table(notas4) #Fr. Abs
prop.table(table(notas4)) #Fr. Rel


## -----------------------------------------------------
cumsum(table(notas4)) #Fr. Abs. Cum
cumsum(prop.table(table(notas4))) #Fr. Rel. Cum


## -----------------------------------------------------
notasHist = hist(notas, breaks = L, right = FALSE, include.lowest = TRUE, 
                 plot = FALSE)
FAbs = notasHist$count
FRel = prop.table(FAbs)
FAbsCum = cumsum(FAbs)
FRelCum = cumsum(FRel)


## -----------------------------------------------------
intervalos = c("[0,5)","[5,7)","[7,9)","[9,10]")
calificacion = c("Suspenso", "Aprobado", "Notable", "Excelente")
marcas = notasHist$mids
tabla.Fr = data.frame(intervalos,calificacion,marcas,FAbs,FAbsCum,
                      FRel,FRelCum)
tabla.Fr


## -----------------------------------------------------
TablaFrecs.L(notas, L, TRUE)


## ---- echo = FALSE------------------------------------
L = c(L1,L2,L3,L4,L5,L6,L7,L8,L9,L10,L11)
L
intervals = as.character(c("[20.95,22.25)","[22.25,23.55)","[23.55,24.85)","[24.85,26.15)","[26.15,27.45)","[27.45,28.75)","[28.75,30.05)","[30.05,31.35)","[31.35,32.65)","[32.65,33.95)"))
TF.L = function(x,L,V){
  x_cut = cut(x, breaks=L, right=FALSE, include.lowest=V)
  mc = (L[1:(length(L)-1)]+L[2:length(L)])/2
  Fr.abs = as.vector(table(x_cut)) 
  Fr.rel = round(Fr.abs/length(x),4)
  Fr.cum.abs = cumsum(Fr.abs)
  Fr.cum.rel = cumsum(Fr.rel)
  tabla = data.frame(intervals, mc, Fr.abs, Fr.cum.abs, Fr.rel, Fr.cum.rel)
  tabla
  }
tabla = TF.L(cw,L,FALSE)
tabla


## -----------------------------------------------------
TOT = tabla$Fr.cum.abs[10]
TOT
anchura.media = round(sum(tabla$Fr.abs*tabla$mc)/TOT,3)
anchura.media #Media
anchura.var = round(sum(tabla$Fr.abs*tabla$mc^2)/TOT-anchura.media^2,3)
anchura.var #Varianza


## -----------------------------------------------------
anchura.dt = round(sqrt(anchura.var),3)
anchura.dt #Desviación típica
I.modal = tabla$intervals[which(tabla$Fr.abs == max(tabla$Fr.abs))]
I.modal #Intervalo modal


## -----------------------------------------------------
I.critic = tabla$intervals[which(tabla$Fr.cum.rel >= 0.5)]
I.critic[1] #Intervalo critic


## -----------------------------------------------------
n = TOT
Lc = L[4]
Lc.pos = L[5]
Ac = L[5]-L[4]
Nc.ant = tabla$Fr.cum.abs[3]
nc = tabla$Fr.abs[4]
M = Lc+Ac*((n/2)-Nc.ant)/nc
M #Aproximación de la mediana de los datos "reales"
median(cw) #Mediana de los datos "reales"


## -----------------------------------------------------
aprox.quantile.p = function(Lcrit,Acrit,n,p,Ncrit.ant,ncrit){
  round(Lcrit+Acrit*(p*n-Ncrit.ant)/ncrit,3)
}
aprox.quantile.p(Lc,Ac,n,0.25,Nc.ant,nc) #Primer cuartil
aprox.quantile.p(Lc,Ac,n,0.75,Nc.ant,nc) #Tercer cuartil


## -----------------------------------------------------
quantile(cw,0.25)
quantile(cw,0.75)


## ---- echo = FALSE------------------------------------
par(mfrow = c(1,2))
hist(cw, breaks = L, right = FALSE, main = "Histograma con intervalos \n de misma anchura",xlab = "")
Lnotas = c(0,5,7,9,10)
hist(notas, breaks = Lnotas, right = FALSE, include.lowest = TRUE, main = "Histograma con intervalos \n de diferente anchura", xlab = "")
par(mfrow = c(1,1))


## -----------------------------------------------------
histAbs = function(x,L) {
  h = hist(x, breaks = L, right = FALSE, freq = FALSE,
           xaxt = "n", yaxt = "n", col = "lightgray", 
           main = "Histograma de frecuencias absolutas", 
           xlab = "Intervalos y marcas de clase",
           ylab = "Frecuencias absolutas")
  axis(1, at=L)
  text(h$mids, h$density/2, labels=h$counts, col="purple") 
  }


## ---- echo=FALSE--------------------------------------
set.seed(1)
edades = c(sample(0:99,80,replace = TRUE),rep(35,10),rep(22,5),rep(17,3),50,50)
extremos = c(0,20,40,60,80,100)
par(mfrow=c(1, 2))
histAbs(edades, extremos)
rug(edades)
histAbs(edades, extremos)
rug(jitter(edades))
par(mfrow=c(1,1))
set.seed(NULL)


## -----------------------------------------------------
histAbsCum = function(x,L) {
  h = hist(x, breaks = L, right = FALSE , plot = FALSE) 
  h$density = cumsum(h$density)
  plot(h, freq = FALSE, xaxt = "n", yaxt = "n", col = "lightgray", 
       main = "Histograma de frecuencias\nabsolutas acumuladas", 
       xlab = "Intervalos", ylab = "Frec. absolutas acumuladas")
  axis(1, at=L)
  text(h$mids, h$density/2, labels = cumsum(h$counts), col = "purple") 
  }


## ---- echo=FALSE--------------------------------------
set.seed(1)
edades = c(sample(0:99,80,replace = TRUE),rep(35,10),rep(22,5),rep(17,3),50,50)
extremos = c(0,20,40,60,80,100)
par(mfrow=c(1, 2))
histAbsCum(edades, extremos)
rug(edades)
histAbsCum(edades, extremos)
rug(jitter(edades))
par(mfrow=c(1,1))
set.seed(NULL)


## ---- echo=FALSE,out.width="50%",fig.align='center'----
knitr::include_graphics("Imgs/gauss.png",dpi=600)


## -----------------------------------------------------
histRel = function(x,L) {
  h = hist(x, breaks=L, right=FALSE , plot=FALSE)
  t = round(1.1*max(max(density(x)[[2]]),h$density),2) 
  plot(h, freq = FALSE, col = "lightgray", 
       main = "Histograma de frec. relativas\ny curva de densidad estimada", 
       xaxt="n", ylim=c(0,t), xlab="Intervalos", ylab="Densidades")
  axis(1, at = L) 
  text(h$mids, h$density/2, labels = round(h$counts/length(x),2), 
       col = "blue")
  lines(density(x), col = "purple", lwd = 2) 
  }


## ---- echo=FALSE--------------------------------------
set.seed(1)
edades = c(sample(0:99,80,replace = TRUE),rep(35,10),rep(22,5),rep(17,3),50,50)
extremos = c(0,20,40,60,80,100)
par(mfrow=c(1, 2))
histRel(edades, extremos)
rug(edades)
histRel(edades, extremos)
rug(jitter(edades))
par(mfrow=c(1,1))
set.seed(NULL)


## -----------------------------------------------------
histRelCum = function(x,L){
  h = hist(x, breaks = L, right = FALSE , plot = FALSE)
  h$density = cumsum(h$counts)/length(x)
  plot(h, freq = FALSE, 
      main = "Histograma de frec. rel. acumuladas\n y 
      curva de distribución estimada", 
      xaxt = "n", col = "lightgray", xlab = "Intervalos", 
      ylab = "Frec. relativas acumuladas") 
  axis(1, at = L)
  text(h$mids, h$density/2, labels = round(h$density ,2), col = "blue")
  dens.x = density(x)
  dens.x$y = cumsum(dens.x$y)*(dens.x$x[2]-dens.x$x[1]) 
  lines(dens.x,col = "purple",lwd = 2)
}


## ---- echo=FALSE--------------------------------------
set.seed(1)
edades = c(sample(0:99,80,replace = TRUE),rep(35,10),rep(22,5),rep(17,3),50,50)
extremos = c(0,20,40,60,80,100)
par(mfrow=c(1, 2))
histRelCum(edades, extremos)
rug(edades)
histRelCum(edades, extremos)
rug(jitter(edades))
par(mfrow=c(1,1))
set.seed(NULL)


## ----eval=FALSE---------------------------------------
## hist(cw, breaks = L, right = FALSE, main = "Histograma de las
##      anchuras de los cangrejos")


## ----echo=FALSE,out.width="50%",fig.align='center'----
hist(cw, breaks = L, right = FALSE, main = "Histograma de las 
     anchuras de los cangrejos")


## -----------------------------------------------------
hist(cw, breaks = L, right = FALSE, plot = FALSE)


## ----out.width="50%",fig.align='center'---------------
histAbs(cw,L)


## ----out.width="50%",fig.align='center'---------------
histAbsCum(cw,L)


## ----out.width="50%",fig.align='center'---------------
histAbs(cw,L)
rug(cw)


## ----out.width="50%",fig.align='center'---------------
histAbs(cw,L)
rug(jitter(cw))


## -----------------------------------------------------
str(density(cw))


## ----out.width="45%",fig.align='center'---------------
histRel(cw,L)


## ---- results="hide", eval = FALSE--------------------
## histRel(cw,L)
## curve(dnorm(x, mean(cw), sd(cw)), col="cyan4", lty=4, lwd=2,
## add=TRUE)
## legend("topright", lwd=c(2,2), lty=c(1,4), col=c("purple","cyan4"),
##        legend=c("densidad estimada","densidad normal"))


## ---- echo=FALSE--------------------------------------
histRel(cw,L)
curve(dnorm(x, mean(cw), sd(cw)), col="cyan4", lty=4, lwd=2,
add=TRUE)
legend("topright", lwd=c(2,2), lty=c(1,4), col=c("purple","cyan4"), legend=c("densidad estimada","densidad normal"))


## ----out.width="50%",fig.align='center'---------------
histRelCum(cw,L)

