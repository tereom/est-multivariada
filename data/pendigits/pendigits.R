library(stringr)

con <- file('pendigits-orig.tra') 
open(con)
results.list <- list()
digitos <- list()
j <- 1
line <- readLines(con, n = 1, warn = FALSE)

#esto soloamente lee el pirmer segmento
while (!is.null(line)) {
 while(str_sub(line,1,4)!='.SEG'){
   line <- readLines(con, n = 1, warn = FALSE)
 }
  print(line)
  digitos[[j]] <- as.integer(str_split_fixed(line, '"', 3)[2])
   while(line!='.PEN_DOWN'){
        line <- readLines(con, n = 1, warn = FALSE)
   }
  
   dat.list <- list()
    i <- 1
  line <- readLines(con, n = 1, warn = FALSE)
   while(line!='.PEN_UP'){
      dat.list[[i]] <- as.integer(str_split_fixed(str_trim(line), " ", 2))
      i <- i + 1
      line <- readLines(con, n = 1, warn = FALSE)
  }
  results.list[[j]] <- Reduce(rbind, dat.list)
  j <- j + 1
  line <- readLines(con, n = 1, warn = FALSE)
  #print(j)
} 
close(con)

plot(results.list[[100]])
digitos[[100]]

convertir.serie <- function(patron){
  d.1 <- rbind(patron[-1,],c(NA,NA)) - patron
  d.2 <- d.1[-nrow(d.1),]
  angulo <- acos(d.2[,1]/(apply(d.2,1, function(x) sqrt(sum(x^2)))))*sign(d.2[,2])
  velocidad <- apply(d.2,1, function(x) sqrt(sum(x^2)))
   #angulo[is.nan(angulo)] <- 0
  dat.out<- data.frame(angulo=angulo[],velocidad = velocidad[])
  dat.out[!is.na(dat.out$angulo), ]
}



trazar <- function(datos){
  angulos <- datos[,1]
  velocidad <- datos[,2]
  longitud <- length(angulos)+1
  puntos <- data.frame(x=rep(NA,longitud), y=rep(NA, longitud))
  puntos[1, ] <- c(0,0)
  for(i in 2:length(angulos)){
    puntos[i,] <- puntos[i-1,] + velocidad[i-1]*c(cos(angulos[i-1]),sin(angulos[i-1]))
  }
  puntos
}

con.ang <- convertir.serie(results.list[[611]])
plot(trazar(con.ang))
plot(results.list[[611]])

## entrenar HMM para distintos dígitos.
filtro.5 <- sapply(digitos, function(x) x==8)
results.5 <- results.list[filtro.5]
plot(results.5[[1]], xlim=c(80,400), ylim=c(80,400))
plot(results.5[[9]], xlim=c(80,400), ylim=c(80,400))
plot(results.5[[19]], xlim=c(80,400), ylim=c(80,400))
plot(results.5[[29]])

plot(trazar(con.ang <- convertir.serie(results.5[[1]])))
plot(trazar(con.ang <- convertir.serie(results.5[[9]])))
plot(trazar(con.ang <- convertir.serie(results.5[[19]])))
plot(trazar(con.ang <- convertir.serie(results.5[[29]])))


library(plyr)
dat.6 <- ldply(1:700, function(i){
  conv.1 <- convertir.serie(results.5[[i]])
  #conv.2 <- conv.1[!is.na(conv.1)]
  data.frame(angulo=conv.1[,1], velocidad=conv.1[,2], serie=i, tiempo=1:nrow(conv.1))
})
dat.6.na <- dat.6
#dat.6.na <- subset(dat.6, !is.na(angulo))
longitudes <- ddply(dat.6, 'serie', summarise, long=max(tiempo))

library(Hmisc)
dat.6.na$angulo.cut <- factor(cut2(dat.6$angulo, g=8, levels.mean=TRUE))
levs <- as.numeric(levels(dat.6.na$angulo.cut))
library(depmixS4)
mod <- depmix(list(angulo.cut~1, velocidad~1), data = dat.6.na, nstates = 16,
    family =  list(multinomial("identity"), gaussian()),
   ntimes = longitudes$long)
if(FALSE){
getpars(mod)
setpars(mod, value=1:npar(mod))
#setpars(mod, getpars(mod, "fixed"))
pars <- c(unlist(getpars(mod)))
pars[17:272] <- 0
pars[sapply(1:16, function(i){0+17*i})] <- 1/2
pars[sapply(1:15, function(i){1+17*i})] <- 1/2
#pars[sapply(1:14, function(i){2+17*i})] <- 1/3
pars[255:256] <- 0.5
pars[272] <- 1
pars[1:16] <- 0
pars[1] <- 1
mod <- setpars(mod, pars)
free <- rep(0,length(pars))
free[1:16] <- 1
free[sapply(1:16, function(i){0+17*i})] <- 0
free[sapply(1:15, function(i){1+17*i})] <- 1
free[sapply(1:14, function(i){2+17*i})] <- 1

free[273:length(pars)] <- 1
}
 

fm <- fit(mod,  emcontrol=em.control(maxit=60))#, fixed=!free)
class(fm) <- 'depmix'



sim.1 <- simulate(fm)
#sim.1@states


resp.1 <- sim.1@response

estados <- sim.1@states
#estados.2 <- sim.1@states[38:(38+63-1)]
angulos.sim <- lapply(resp.1[estados], function(x){ 
  params <- x[[1]]@parameters$coefficients
  c(sample(levs, 1, prob=params ), rnorm(1, x[[2]]@parameters[1][[1]], sd=x[[2]]@parameters[2]$sd))
   # rnorm(1, params[[1]], sd=params$sd )}
}
)
dat.x <- data.frame(Reduce(rbind, angulos.sim))
names(dat.x) <- c('angulo','velocidad')
#estados
head(dat.x, 6)
plot(trazar(dat.x[1:77,]), type='l')
plot(trazar(dat.x[78:(77+37),]), type='l')
plot(trazar(dat.x[(44+43):(44+42+43),]), type='l')
plot(trazar(dat.x[1:1000,]), type='l')

##### Con HMM
#library(HMM)
#hmm <- initHMM(letters[1:16], levs, transProbs=matrix(1/16,16,16), emissionProbs=matrix(1/16,16,16))    
#print(hmm)
#fit.bm <- baumWelch(hmm, dat.6.na$angulo.cut, maxIterations=10, delta=1E-9, pseudoCount=0)
#fit.viterbi <- viterbiTraining(hmm, dat.6.na$angulo.cut, maxIterations=20)

