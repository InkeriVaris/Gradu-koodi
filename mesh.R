# Verkon sovitus INLA-paketilla

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
library(INLA)

# Ladataan aineisto
load("fishreef_full.Rdata")

# Tehdään mesh/verkko
dim(dcoords)
summary(dcoords)

# Sovitetaan verkko, eli diskretisoidaan spatiaalinen kenttä
bnd <- inla.nonconvex.hull(dcoords, convex = 0.1, resolution = c(63, 83))
mesh <- inla.mesh.2d(
  loc = dcoords,
  boundary = bnd,
  max.edge = c(0.2, 0.7), 
  cutoff = 0.1 
)

# Kuva verkosta
plot(mesh)

# Matern-SPDE-malli objekti
fishreef_mesh = inla.spde2.matern(mesh, alpha = 2)

# Tallennetaan verkko
save(fishreef_mesh, file = "fishreef_mesh.Rdata")



