install.packages("ncdf4")
library(ncdf4)
install.packages("plotrix")
library(plotrix)

fn = "C:/Users/tryst/Desktop/PhD/Chl_a_Baltic/CHLa_SW_BalticSea_2017_07_01-2017_08_31_Copernicusdatabase.nc"
nc = nc_open(fn)

print(nc)
attributes(nc$var)

lon = ncvar_get(nc, "lon")
nlon = dim(lon)
head(lon)

lat= ncvar_get(nc, "lat", verbose = FALSE)
nlat = dim(lat)
head(lat)
print(lat)
print(lon)

print(c(nlon, nlat))

t = ncvar_get(nc, "time")
print(t)
tunits = ncatt_get(nc, "time", "units")
nt = dim(t)
nt
tunits

CHL_array = ncvar_get(nc, "CHL")
dlname = ncatt_get(nc, "CHL", "long_name")
dlname <- ncatt_get(nc,"CHL","long_name")
dunits <- ncatt_get(nc,"CHL","units")
fillvalue <- ncatt_get(nc,"CHL","_FillValue")
dim(CHL_array)
print(CHL_array)

dim(CHL_array)

####### KIEL
kiel = CHL_array[54:58, 127:131, 1:62]
x = sum(!is.na(kiel))
x
se1 = sd(kiel, na.rm = TRUE)/sqrt(x)
mean(kiel, na.rm = TRUE)
sd(kiel, na.rm = TRUE)
se1

####### AHP
ahp = CHL_array[173:177, 131:135, 1:62]
y = sum(!is.na(ahp))
y
se2 = sd(ahp, na.rm = TRUE)/sqrt(y)
mean(ahp, na.rm = TRUE)
sd(ahp, na.rm = TRUE)
se2

#######USE
use = CHL_array[272:276, 160:164, 1:62]
z = sum(!is.na(ahp))
z
se3 = sd(ahp, na.rm = TRUE)/sqrt(z)
mean(use, na.rm = TRUE)
sd(use, na.rm = TRUE)
se3

print(dat)
dim(dat)
dat[47, 20, 10]
mean(dat, na.rm = TRUE)
dat2 = apply(dat, c(1, 2), mean, na.rm = TRUE)
print(dat2)
mean(dat2, na.rm = TRUE)

######################################## example

Daily_b06_45
fn1 = "C:/Users/tryst/Desktop/PhD/Chl_a_Baltic/Daily_b06_45.nc"
nc1 = nc_open(fn1)
print(nc1)
