# Tämä tiedosto sisältää mallin sovituksen TMB-paketilla
# ja koodin kuvien piirtoa varten

library(TMB)
library(Matrix)
library(INLA)

# Käännetään cpp-tiedosto, tämä luo workin directoryyn .dll -päätteisen tiedoston
compile(
  "latenttimalli.cpp", 
  framework = "TMBad",
  openmp = TRUE
)

# Ladataan edellisen funktion luoma käännetty dll-tiedosto R:ään
dyn.load(dynlib("latenttimalli"))
openmp(max = TRUE, DLL = "latenttimalli", autopar = TRUE)

# Ladataan aineisto
load("fishreef_full.Rdata")

# Ladataan INLA-paketilla luotu verkko-tiedosto 
load("fishreef_mesh.Rdata")

# Kuva verkosta
plot(fishreef_mesh$mesh)

# Luodaan A matriisi, joka yhdistää verkon solmut ja havainnot
A = inla.spde.make.A(
  fishreef_mesh,
  loc = cbind(Xreef$Start_Longitude, Xreef$Start_Latitude)
)

# Verkon solmujen määrä
N = fishreef_mesh$n.spde 

# Standardoidaan jatkuvaluonteiset kovariaatit
Xreef$Start_Depth <- scale(Xreef$Start_Depth)
Xreef$Substrate <- scale(Xreef$Substrate)
Xreef$LastOfTemp <- scale(Xreef$LastOfTemp)

# Kategorisoidaan diskreetit kovariaatit
Xreef$Station_ID <- factor(Xreef$Station_ID)
Xreef$Turbidity <- factor(Xreef$Turbidity)
Xreef$Current_Direction <- factor(Xreef$Current_Direction)


# Alustetaan data ja parametrit
n = nrow(yreef)
p = ncol(yreef)
d = 2 # latenttien muuttujien lukumäärä
Y = as.matrix(yreef)

# Muodostetaan design-matriisi X, joka sisältää ympäristökovariaatit:
X = model.matrix(~ 1 + Start_Depth + Substrate + I(Year -2011) + LastOfTemp + Turbidity + Current_Direction, data = Xreef)

# Lista data-objekteista: Tässä listassa kaikki DATA_-alkuiset objektit, jotka alustetaan .cpp tiedostossa
projmat <- A
data.list <- list(y=Y, 
                  x=X, 
                  spde = fishreef_mesh$param.inla[c("M0", "M1", "M2")],
                  projmat=projmat)

# Lista parametreista: Tässä listassa kaikki PARAMETER-alkuiset objektit, jotka alustetaan .cpp tiedostossa
# Nämä ovat estimoitavia parametreja ja satunnaisvaikutuksia:
parameter.list <- list(
  b = matrix(0, ncol(X), p), #kovariaatit
  u = array(0, dim = c(ncol(projmat), d)), #latentit muuttujat
  lg_phi = rep(log(0.5), times = p), #NB-jakauman hajontaparamteri
  lataus = rep(1, times = p*d -(d*(d-1))/2), #latausmatriisi 
  lg_range = log(0.3), #matern-kovariassin parametri
  lv_cor = rep(0, times = (d*(d-1))/2) #latenttien muuttujien kovarianssi 
  )

# Alustetaan malli: 
# Funktio "MakeADFun" muodostaa Laplace-approksimoidun uskottavuusfunktion ja luo gradienttifunktiot
# Sille annetaan edellä luodut listat data -objekteista ja parametreista:
obj = MakeADFun(
  data = data.list, # lista data-objekteista
  silent = FALSE, 
  parameters = parameter.list, # lista parametreista
  DLL = "latenttimalli", # dll tiedoston nimi
  random = c("u"), # parametrien nimet, jotka ovat satunnais-efektejä, ja joiden suhteen marginaali uskottavuus lasketaan ja Laplace-approksimoidaan
  smartsearch = TRUE
)

# Parametrivektori
obj$par
# Laplace-approksimoitu marginaali-uskottavuus
obj$fn(obj$par) 
# Gradienttifunktio
obj$gr(obj$par)
# Latausmatriisi
obj$report()

# Mallin estimointi gradienttipohjaisella optimointifunktiolla
t1 <- system.time({
  opt = nlminb(
    obj$par,
    obj$fn,
    obj$gr,
    control = list(eval.max = 1000, iter.max = 1000, rel.tol = 1e-10)
  )
})[3]

# Suurimman uskottavuuden estimaatti
-opt$objective
# Estimoidut parametrit
opt$par

# Tarkistetaan konvergoiko optimointi prosessi gradienttien avulla
opt$convergence
# Gradientit pääosin lähellä nollaa, mutta muutamien parametrien osalta ei niin hyvä:
plot(c(obj$gr(opt$par)))

# Estimoidut parametrit, satunnaisvaikutusket ja matriisit
obj$env$parList()

# Parametrien keskivirheet:
sds<- sdreport(obj)
summary(sds)


library(corrplot)
library(gclus)

# Piirretään jäännöskovarianssimatriisi huomioiden ylihajonta NB-jakauman tapauksessa
# Jäännöskovarianssi
ResCov <- t(obj$report()$lam) %*% obj$report()$lam
ResCov <- ResCov + diag(log(exp(obj$env$parList()$lg_phi)+1),
                        ncol = ncol(ResCov))

# Tehdään tästä korrelaatiomatriisi
cormatrix <- cov2cor(ResCov)
colnames(cormatrix) <- colnames(Y)
rownames(cormatrix) <- colnames(Y)

# Piirretään jäännöskovarianssimatriisi
corrplot(cormatrix[order.single(cormatrix), order.single(cormatrix)], diag = F, type = "lower", 
         method = "square", tl.cex = 0.5, tl.srt = 45, tl.col = "black")


library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(Matrix)

# Alustetaan kartta tutkimusalueesta kuvia varten
world <- ne_countries(scale = "medium", returnclass = "sf")

# Koordinaatit
lon_min <- min(dcoords[,1])-3
lon_max <- max(dcoords[,1])+2
lat_min <- min(dcoords[,2])-3
lat_max <- max(dcoords[,2])+1

# Kartta tutkimusalueesta
map <- ggplot(data = world) +
  geom_sf() +
  coord_sf(xlim = c(lon_min, lon_max), ylim = c(lat_min, lat_max), expand = FALSE) +
  theme_minimal()
map

# Piirretään kalalajien spatiaaliset jakaumat tutkimusalueen kartalle

# Yhdistetään verkon solmut koordinaatteihin
Ad = inla.spde.make.A(fishreef_mesh, dcoords)

# Latentit muuttujat kertaa lajien latausmatriisi
LVgamma <- Ad %*% obj$env$parList()$u %*% obj$report()$lam

# Koordinaatit ja latentit muuttujat kertaa lajien latausmatriisi
dcoords2 <- data.frame(V1 = dcoords[,1], V2 = dcoords[,2], LVgamma = as.matrix(LVgamma))


# Kaikki lajien spatiaaliset jakaumat tutkimusalueella samassa kuvassa
lajit <- data.frame(ryhma = c(rep("Balistes capriscus", times = nrow(dcoords2)), rep("Caulolatilus microps", times = nrow(dcoords2)), rep("Centropristis striata", times = nrow(dcoords2)),rep("Cephalopholis cruentata", times = nrow(dcoords2)),
                              rep("Epinephelus adscensionis", times = nrow(dcoords2)),rep("Epinephelus itajara", times = nrow(dcoords2)), rep("Epinephelus morio", times = nrow(dcoords2)),rep("Epinephelus niveatus", times = nrow(dcoords2)),
                              rep("Haemulon plumierii", times = nrow(dcoords2)),rep("Lachnolaimus maximus", times = nrow(dcoords2)),rep("Lutjanus analis", times = nrow(dcoords2)),rep("Lutjanus campechanus", times = nrow(dcoords2)),rep("Lutjanus griseus", times = nrow(dcoords2)),
                              rep("Malacanthus plumieri", times = nrow(dcoords2)),rep("Mycteroperca interstitialis", times = nrow(dcoords2)),rep("Mycteroperca microlepis", times = nrow(dcoords2)),rep("Mycteroperca phenax", times = nrow(dcoords2)),
                              rep("Ocyurus chrysurus", times = nrow(dcoords2)),rep("Pagrus pagrus", times = nrow(dcoords2)),rep("Rhomboplites aurorubens", times = nrow(dcoords2)),rep("Seriola dumerili", times = nrow(dcoords2))) 
                       ,LVgamma = c(dcoords2$LVgamma.1, dcoords2$LVgamma.2,dcoords2$LVgamma.3,dcoords2$LVgamma.4, dcoords2$LVgamma.5, dcoords2$LVgamma.6, dcoords2$LVgamma.7, dcoords2$LVgamma.8, dcoords2$LVgamma.9, dcoords2$LVgamma.10, dcoords2$LVgamma.11,
                                    dcoords2$LVgamma.12,dcoords2$LVgamma.13, dcoords2$LVgamma.14, dcoords2$LVgamma.15, dcoords2$LVgamma.16, dcoords2$LVgamma.17, dcoords2$LVgamma.18, dcoords2$LVgamma.19, dcoords2$LVgamma.20, dcoords2$LVgamma.21),
                       V1 = c(rep(dcoords2$V1, times = 21)),
                       V2 = c(rep(dcoords2$V2, times = 21)))

map +
  geom_point(data = lajit, aes(x = V1, y = V2, colour=LVgamma), size = 1) +
  labs(title = "Kalalajit", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LM") + 
  facet_wrap(~ryhma)

# Lajien Centropristis striata, Balistes capriscus ja Mycteroperca phenax spatiaaliset jakaumat tutkimusalueella samassa kuvassa
osalajit <- data.frame(ryhma = c(rep("Centropristis striata", times = nrow(dcoords2)),
                              rep("Balistes capriscus", times = nrow(dcoords2)),
                              rep("Mycteroperca phenax", times = nrow(dcoords2))) 
                    ,LVgamma = c(dcoords2$LVgamma.3, dcoords2$LVgamma.1,dcoords2$LVgamma.17),
                    V1 = c(rep(dcoords2$V1, times = 3)),
                    V2 = c(rep(dcoords2$V2, times = 3)))

map +
  geom_point(data = osalajit, aes(x = V1, y = V2, colour=LVgamma), size = 1) +
  labs(x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LM*lataukset") + 
  facet_wrap(~ryhma)


# Yksittäiset kuvat jokaisen lajin spatiaalisesta jakaumasta tutkimusalueella
map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.1), size = 2) +
  labs(title = "Balistes capriscus", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings") 

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.2), size = 2) +
  labs(title = "Caulolatilus microps", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.3), size = 2) +
  labs(title = "Centropristis striata", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.4), size = 2) +
  labs(title = "Cephalopholis cruentata", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.5), size = 2) +
  labs(title = "Epinephelus adscensionis", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.6), size = 2) +
  labs(title = "Epinephelus itajara", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.7), size = 2) +
  labs(title = "Epinephelus morio", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.8), size = 2) +
  labs(title = "Epinephelus niveatus", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.9), size = 2) +
  labs(title = "Haemulon plumierii", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.10), size = 2) +
  labs(title = "Lachnolaimus maximus", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.11), size = 2) +
  labs(title = "Lutjanus analis", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.12), size = 2) +
  labs(title = "Lutjanus campechanus", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.13), size = 2) +
  labs(title = "Lutjanus griseus", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.14), size = 2) +
  labs(title = "Malacanthus plumieri", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.15), size = 2) +
  labs(title = "Mycteroperca interstitialis", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.16), size = 2) +
  labs(title = "Mycteroperca microlepis", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.17), size = 2) +
  labs(title = "Mycteroperca phenax", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.18), size = 2) +
  labs(title = "Ocyurus chrysurus", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.19), size = 2) +
  labs(title = "Pagrus pagrus", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.20), size = 2) +
  labs(title = "Rhomboplites aurorubens", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")

map +
  geom_point(data = dcoords2, aes(x = V1, y = V2, colour=LVgamma.21), size = 2) +
  labs(title = "Seriola dumerili", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV*loadings")


# Piirretään kuvat molemmille latenteille muuttujille

# Latentit muuttujat
LVgammaLat <- Ad %*% obj$env$parList()$u 

# Koordinaatit ja latentit muuttujat
dcoords3 <- data.frame(V1 = dcoords[,1], V2 = dcoords[,2], LVgamma = as.matrix(LVgammaLat))

# Latentin muuttujan 1 kuva
map +
  geom_point(data = dcoords3, aes(x = V1, y = V2, colour=LVgamma.1), size = 2) +
  labs(title = "Latentti muuttuja 1", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV")

# Latentin muuttujan 2 kuva
map +
  geom_point(data = dcoords3, aes(x = V1, y = V2, colour=LVgamma.2), size = 2) +
  labs(title = "Latentti muuttuja 2", x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LV")


# Molemmat latentit muuttujat samassa kuvassa
latentit <- data.frame(ryhma = c(rep("Latentti muuttuja 1", times = nrow(dcoords3)), rep("Latentti muuttuja 2", times = nrow(dcoords3))), 
                       LVgamma = c(dcoords3$LVgamma.1, dcoords3$LVgamma.2),
                       V1 = c(dcoords3$V1, dcoords3$V1),
                       V2 = c(dcoords3$V2, dcoords3$V2))

map +
  geom_point(data = latentit, aes(x = V1, y = V2, colour=LVgamma), size = 2) +
  labs(x = "Pituusaste", y = "Leveysaste") + scale_colour_gradientn(colours = hcl.colors(5, palette = "inferno"), name = "LM") + 
  facet_wrap(~ryhma)


library("lattice")

# Piirretään yhteinen kuva kovariaateista

# Parametrien keskivirheet
par_sd <- (summary(sds))

# Ympäristökovariaatien estimaatit
coeftab <- obj$env$parList()$b
colnames(coeftab) = colnames(obj$env$data$y)
rownames(coeftab) = colnames(obj$env$data$x)

# Estimaattien keskihajonta
coefSd <- matrix(par_sd[rownames(par_sd) == "b",2], nrow(coeftab), p)

# 95%:n luottamusvälit
CIlow <- coeftab + qnorm(0.025)*coefSd
CIupp <- coeftab + qnorm(0.975)*coefSd
Asign <- ((CIlow>0) | (CIupp<0))*sign(coeftab)

# Asetetaan ympäristökovariaattien nimet
varnames <- c("Vakiotermi","Mittauspisteen syvyys","Merenpohjan alusta","Vuosi","Lämpötila","Melko selkeä vesi","Samea vesi","Virtaus sivulle","Virtaus kohti")

expT =FALSE

# Funktio kuvan piirtoa varten
{
  rownames(coeftab) <- varnames
  # Ei piirretä vakiotermiä:
  Asign1 = Asign[-1,]
  if(expT) {
    Aplot <- exp(coeftab[-1,])
  } else {
    Aplot = coeftab[-1,]
  }
  y <- 1:ncol(Aplot)
  x <- 1:nrow(Aplot)
  
  grid <- expand.grid(x=x, y=y)
  grid$z <- c((Aplot[,ncol(Aplot):1]))
  w <- grid$w <- abs(c((Asign1)))
  a <- max(abs(c(Aplot)))
  
  # Väriasteikko:
  # Väriasteikon keskusta
  keskus = 0
  if(expT){
    colort <- colorRampPalette(c("royalblue1", "lightskyblue", "white", "red", "red3"))
    p <- levelplot((Aplot), grid, 
                   panel=function(...) {
                     arg <- list(...)
                     panel.levelplot(...)
                     panel.text(arg$x, arg$y, ifelse(w>0,round(arg$z,2),""))
                   }, 
                   at = exp(c(seq(keskus - log(a), keskus , length = 25)[-25] ,seq(keskus, (keskus + log(a)), length = 25))), ylab = "", main="", scales=list(x=list(rot=45, cex=0.8), y=list(cex=0.8)),
                   xlab = "", col.regions = colort(100)) 
    
  } else {
    colort <- colorRampPalette(c("royalblue1", "lightskyblue", "white", "red", "red3"))
    p <- levelplot((Aplot), grid, 
                   at = c(c(-a,-15,-10,seq((keskus - 5 ), (keskus + 5), length = 51),10,15,a)), ylab = "", main="", scales=list(x=list(rot=45, cex=0.8), y=list(cex=0.8)),
                   xlab = "", col.regions = colort(100)) 
  }
  
  print(p)
}

dev.off()


# Piirretään kovariaattien estimaatit ja 95% luottamusvälit jokaiselle kalalajille

# Ympäristökovariaattien estimaatit
coeftab <- obj$env$parList()$b 
colnames(coeftab) = colnames(obj$env$data$y)
rownames(coeftab) = colnames(obj$env$data$x)
Xcoef <- t(coeftab)

# Estimaattien keskihajonnat
kertoimet <- matrix(sqrt(diag(sds$cov.fixed))[1:189],9,21)
sdXcoef  <- t(kertoimet)

# Asetetaan kovariaattien nimet
cnames <- c("Vakiotermi","Mittauspisteen syvyys","Merenpohjan alusta","Vuosi","Lämpötila","Melko selkeä vesi","Samea vesi","Virtaus sivulle","Virtaus kohti")
labely <- rownames(Xcoef)
k <- length(cnames)
m <- length(labely)

# Asetetaan x-akselille tarvittavat rajat kovariaatille "Samea vesi"
xlim.list <- list(NULL,NULL,NULL,NULL,NULL,NULL,c(-50,50),NULL,NULL)

# Kovariaattien estimaatit ja 95%:n luottamusvälit jokaiselle lajille
# Kolme kuvaa samalla sivulla:
par(mfrow = c(1, 3), mai=c(0.8,1.2,0.6,0.1))
for (i in 1:k) {
  Xc <- Xcoef[, i]
  if(nrow(Xcoef)<2) names(Xc) <- rownames(Xcoef)
  lower <- Xc - 1.96 * sdXcoef[, i]
  upper <- Xc + 1.96 * sdXcoef[, i]
  lower <- lower[names(Xc)]
  upper <- upper[names(Xc)]
  
  col.seq <- rep("black", m)
  col.seq[lower < 0 & upper > 0] <- "grey"
  
  At.y <- seq(1, m)
   if (!is.null(xlim.list[[i]])) {
     plot( x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = xlim.list[[i]], pch = "x", cex.lab = 1 )
   } else {
  
    plot( x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = c(min(lower), max(upper)), pch = "x", cex.lab = 1 )
  }
  
  segments( x0 = lower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq )
  abline(v = 0, lty = 1)
  axis( 2, at = At.y, labels = names(Xc), las = 1, cex.axis = 0.6)
}

# Kovariaattien estimaatit ja 95%:n luottamusvälit jokaiselle lajille
# Kaksi kuvaa samalla sivulla ilman vakiotermiä:
par(mfrow = c(1, 2), mai=c(0.6,1.4,0.6,0.1))
for (i in 2:k) {
  Xc <- Xcoef[, i]
  if(nrow(Xcoef)<2) names(Xc) <- rownames(Xcoef)
  lower <- Xc - 1.96 * sdXcoef[, i]
  upper <- Xc + 1.96 * sdXcoef[, i]
  lower <- lower[names(Xc)]
  upper <- upper[names(Xc)]
  
  col.seq <- rep("black", m)
  col.seq[lower < 0 & upper > 0] <- "grey"
  
  At.y <- seq(1, m)
  if (!is.null(xlim.list[[i]])) {
    plot( x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = xlim.list[[i]], pch = "x", cex.lab = 1 )
  } else {
    
    plot( x = Xc, y = At.y, yaxt = "n", ylab = "", col = col.seq, xlab = cnames[i], xlim = c(min(lower), max(upper)), pch = "x", cex.lab = 1 )
  }
  
  segments( x0 = lower, y0 = At.y, x1 = upper, y1 = At.y, col = col.seq )
  abline(v = 0, lty = 1)
  axis( 2, at = At.y, labels = names(Xc), las = 1, cex.axis = 0.6)
}

