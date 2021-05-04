install.packages("oceanmap")
install.packages("ncdf4")
install.packages("raster")
install.packages("viridis")

library(oceanmap)
library(ncdf4)
library(raster)
library(viridis)

fn = "dataset-bal-analysis-forecast-phy-monthlymeans_1600248953439.nc"
nc = nc_open(fn)

n <- colorRampPalette(c("blue", "aquamarine4", "khaki1"))(100)


v(fn, pal=n, cbpos ="r")
v ( nfiles [fn], cbpos ="r", replace . na= TRUE )

data("cmap") # load color maps data
names(cmap) # list available color maps

for ( n in names ( cmap ) ) v ( gz . files [2], v _ area =lion , subplot =TRUE ,
pal =n , adaptive . vals =TRUE , main =n )

print(nc)
attributes(nc$var)
