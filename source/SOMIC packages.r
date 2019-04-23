my_packages = c('data.table', 'plyr', 'tools', 'Rcpp', 'lubridate',
                'ggplot2', 'cowplot', 'tcltk2', 'scales',
                'maps', 'maptools', 'rgdal', 'raster', 'rasterVis', 'ncdf4',  'gdalUtils',
                'SoilR', 'zoo','regSmooth', 'graticule')

for (package in my_packages)  {
  if (!require(package, character.only=T, quietly=T)) {
    install.packages(package)
    library(package, character.only=T)
  }
}

rasterOptions(maxmemory = 1e+08, progress = "text")
