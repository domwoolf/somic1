#default values for model parameters
fit.par = c(mic_vmax = 1.1086, mic_km = 0.4089, kdissolution = 2.4893, kdepoly = 0.0436, 
            kdeath_and_exudates = 0.2357, kdesorb = 0.0055, ksorb = 0.3066, 
            kmicrobial_uptake = 1.4506, cue_0 = 0.2811, mcue = 0.0081, mclay = 0.0145)

fit.lower = c(mic_vmax = 1.0, mic_km = 0.1, kdissolution = 0.3, kdepoly = 0.008,
              kdeath_and_exudates = 0.02, kdesorb = 0.0005, ksorb = 0.1, kmicrobial_uptake = 0.1,
              cue_0 = 0.15, mcue = 0, mclay = 0.0)
fit.upper = c(mic_vmax = 2.0, mic_km = 0.6, kdissolution = 6, kdepoly = 1.5,
              kdeath_and_exudates = 1.0, kdesorb = 0.01, ksorb = 1.0, kmicrobial_uptake = 2.0,
              cue_0 = 0.6, mcue = 0.02, mclay = 1/23)

training.set = c("bb3", "expt_pendleton_nB_N90", "expt_pendleton_sB_N45", "expt_sanborn_cwn", 
  "expt_pendleton_nB_N45", "expt_sanborn_cwm", "hoos7_2", "expt_pendleton_sB_N0", 
  "expt_pendleton_nB_MN", "expt_pendleton_sB_N90", "expt_sanborn_cmn")

if (time_step == 'day') {
  fit.par[3:8] = fit.par[3:8] * 12/365
  fit.lower[3:8] = fit.lower[3:8] * 12/365
  fit.upper[3:8] = fit.upper[3:8] * 12/365
}

par.day.to.month = function(){
    fit.par[3:8] = fit.par[3:8] / (12/365)
    fit.lower[3:8] = fit.lower[3:8] / (12/365)
    fit.upper[3:8] = fit.upper[3:8] / (12/365)
}

# function to evaluate list of expressions on data.table, sequentially
with.dt = function(dt, expr){
  for (j in 1:length(expr)) set(dt, NULL, names(expr)[j], dt[, eval(expr[[j]])])
}

# function to calculate temperature dependence of decomposition rate
fT.PrimC = function(temp, method = 'Century2', t.ref = 30){
  switch(method,
    RothC = fT.RothC(temp),
    Century1 = fT.Century1(temp) * fT.RothC(t.ref) / fT.Century1(t.ref),
    Century2 = fT.Century2(temp) * fT.RothC(t.ref) / fT.Century2(t.ref),
    stop('Unrecognised temperature function'))
}

wsat = function(sand, silt, clay) 0.6658*silt + 0.1567*sand - 0.0079*silt^2 - 12.31121/sand -
  6.4756*log(sand) - 0.0038*clay*silt + 0.0038*clay*sand - 0.0042*silt*sand + 52.7526

wfield = function(sand, silt, clay) 118.932*clay + 119.0866*silt + 119.1104*sand + 162.31731/clay -
  46.21921/silt-5.12991/sand  + 18.1733*log(clay) + 0.0013*clay*silt + 0.0022*silt*sand - 11939.3493

wwilt = function(sand, silt, clay) 92.3851 - 1.5722*silt - 0.5423*sand - 0.0072*clay^2 + 0.0072*silt^2 -
  0.0059*sand^2 + 160.14591/clay  +  6.60011/sand + 0.0022*clay*silt - 0.0039*clay*sand

unrotate = function(x) { # inverse of rotate: convert to 0 to 360 longitudes
  raster::shift(rotate(raster::shift(x, 180)), 180)
}

pl = function(x){
  xname = deparse(substitute(x))
  plot(rotate(x), ylim=c(-90,90), main = xname)
  lines(world, lwd=0.3)
}

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

sf = scale_fill_gradientn(
  colours = myPalette(50),
  trans = "log",
  breaks = (c(1,2,5,10,20,40,80)),
  limits = (c(0.99,80)),
  na.value = NA,
  name = expression("SOC (kg m"^-2*")"))

t <- theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position=c(0.42,-0.1),
          legend.title.align=1,
          legend.title = element_text(size=12),
          legend.text = element_text(size = 12),
          legend.direction="horizontal",
          legend.justification = c(0.4,0),
          plot.margin = unit(c(0,0,0,0), "cm"))

g.map = function (r, map_theme=t, map_cols=sf){
  r[r<1] = 1
  gplot((r)) +
    geom_path(data = grat, mapping=aes(long,lat, group=group), size=0.3, linetype=2, color='grey50') +
    geom_raster(aes(fill=value)) +
    geom_path(data = wmap_df, mapping = aes(long, lat, group=group), size=0.3) +
    geom_path(data=bbox_df, mapping=aes(long,lat, group=group), size=0.3) +
    guides(fill= guide_colorbar(barwidth=12, title.vjust=1, nbin=500, draw.ulim = FALSE, draw.llim = FALSE)) +
    coord_equal() +
    map_cols +
    map_theme
}

elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

get_data = function(data_path) {
  files <- list.files(path=data_path, pattern="*.socdat", full.names=T, recursive=FALSE)
  experiments <<- basename(file_path_sans_ext(files))
  all.data = rbindlist(fill = TRUE,
    lapply(files, function(fn) 
      fread(fn, colClasses = c(added.mb='numeric'))[, 
        exp := basename(file_path_sans_ext(fn))]))
  all.data <- merge (all.data, exp.const, by="exp")
  setnames(all.data, "init.soc", "soc")
  all.data
}


Mode <- function(x, na.rm = FALSE) {
  if (na.rm) x = x[!is.na(x)]
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

initialise.somic.data <- function(soc.data) {
  # function to intialise values in soc dataframe
  with.dt(soc.data, expression(
    spm = soc * 0.001,
    ipm = soc * 0.081,
    doc = 0.0,
    mb  = soc * 0.0331,
    mac = soc * (1-0.001-0.081-0.0331),
    atsmd = 0.0,
    sat = 0.0,
    h2o = 0.0,
    a = fT.PrimC(temp, method = 'Century1', t.ref = 28),
    c = ifelse (cover==0, 1, 0.6),
    mic = 1,
    added.doc = 0.0,
    velocity = 0.0,
    add_14c_age = 0.0,
    spm.d13c = soc.d13c,
    ipm.d13c = soc.d13c,
    doc.d13c = soc.d13c,
    mb.d13c  = soc.d13c,
    mac.d13c = soc.d13c,
    co2.d13c = 0.0
  ))
  soc.data
}

# function to intialise values in soc dataframe
initialise.daily.data <- function(soc.data) {
  #convert monthly to daily data
  start.date = as.Date('1800-01-01')
  soc.data[, date := start.date %m+% months(time-1)]
  daily.soc.data = soc.data[, .(date = seq(date[1], date[.N], by='1 day')), by = exp]
  daily.soc.data = merge(soc.data, daily.soc.data, by=c('exp', 'date'), all.y = T)
  setorder(daily.soc.data, exp, date)
  
  #interpolate missing daily values
  daily.soc.data[, time := 1:.N, by = exp]
  daily.soc.data[, date := NULL]
  daily.soc.data[, temp := na.spline(temp), by = exp]
  daily.soc.data[, precip := na.spline(precip)*12/365, by = exp]
  daily.soc.data[, pet := na.spline(pet)*12/365, by = exp]
  daily.soc.data[, cover := na.locf(cover), by = exp]
  daily.soc.data[, clay := na.locf(clay), by = exp]
  daily.soc.data[, depth := na.locf(depth), by = exp]
  daily.soc.data[, soc := na.locf(soc), by = exp]
  daily.soc.data[, soc.d13c := na.locf(soc.d13c), by = exp]
  daily.soc.data[, max_tsmd := na.locf(max_tsmd), by = exp]
  daily.soc.data[is.na(added.spm), added.spm := 0.0]
  daily.soc.data[is.na(added.ipm), added.ipm := 0.0]
  daily.soc.data[is.na(added.mb), added.mb := 0.0]
  daily.soc.data[is.na(added.mac), added.mac := 0.0]
  daily.soc.data[, added.d13c := na.locf(added.d13c), by = exp]
  daily.soc.data = initialise.somic.data(daily.soc.data)
  daily.soc.data
}

initialise.rothc.data <- function(soc.data) {
  iom = soc.data[1, 0.049 * soc^1.139]
  soc.data[, soc := soc - iom]
  soc.data[, dpm := soc * 0.001]
  soc.data[, rpm := soc * 0.081]
  soc.data[, bio := soc * 0.0331]
  soc.data[, sta := soc * (1-0.001-0.081-0.0331)]
  soc.data[, atsmd := 0.0]
  soc.data[, a := fT.RothC(temp)]
  soc.data[, c := ifelse (cover==0, 1, 0.6)]
  soc.data
}

wmean.bricks = function(b.list, w) {
  b.wmean = b.list[[1]] * w[1]
  for (i in seq_along(w)[-1]){
    b.wmean = b.wmean + b.list[[i]] * w[i]
  }
  b.wmean / sum(w)
}
