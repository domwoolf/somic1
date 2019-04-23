# This file demonstrates use of the SOMic model 
# using data from long-term experimental stations

#############################################
#       Initialisation
#############################################
library(rstudioapi)
base_path = dirname(dirname(getActiveDocumentContext()$path))
R_path = paste0(base_path, '/Layers')
data_path = paste0(base_path, '/data')
results_path = paste0(base_path, '/results/')

time_step = 'month'
train = FALSE
save_images = FALSE

setwd(R_path)
package_file = sort(list.files(pattern = 'packages'), decreasing = TRUE)[1]
functions_file = sort(list.files(pattern = 'Functions'), decreasing = TRUE)[1]
cpp_file = sort(list.files(pattern = 'SOMIC.+cpp'), decreasing = TRUE)[1]
source(package_file)
source(functions_file)
sourceCpp(cpp_file)

exp.const = fread(paste0(data_path, "/experiment_constants.csv"))
exp.const[, max_tsmd := -(20+1.3*clay - 0.01*clay^2) * depth/23]

raw.soc.data = get_data(data_path)
all.data = initialise.somic.data(copy(raw.soc.data))

#############################################
#      Solve for optimal parameter values
#############################################
if (train){
  scam_objective = function(p.optim, p.complete, dt, experiments) {
    p = p.complete
    p[names(p.optim)] = p.optim
    use_atsmd = TRUE
    monthly = TRUE
    for (experiment in experiments) {
      max_tsmd = exp.const[exp==experiment, max_tsmd]
      clay = exp.const[exp==experiment, clay]
      somic.out = as.data.table(somic(dt[exp==experiment],
        mic_vmax = p['mic_vmax'],
        mic_km = p['mic_km'],
        kdissolution = p['kdissolution'],
        kdepoly = p['kdepoly'],
        kdeath_and_exudates = p['kdeath_and_exudates'],
        kdesorb = p['kdesorb'],
        ksorb = p['ksorb'],
        kmicrobial_uptake = p['kmicrobial_uptake'],
        cue_0 = p['cue_0'],
        mcue = p['mcue'],
        mclay = p['mclay'],
        clay = exp.const[exp==experiment, clay], 
        max_tsmd = exp.const[exp==experiment, max_tsmd],
        use_atsmd = TRUE))
      dt[exp==experiment, soc := somic.out$soc]
    }
    return(dt[, sum((soc - measured.soc)^2, na.rm = TRUE)])
  }

  training.set = sample(experiments, length(experiments)%/% 2)
  train.data = all.data[exp %in% training.set]
  p = fit.par[!(names(fit.par) %in% c(''))]
  fit = optim(p, scam_objective, p.complete = fit.par, dt = train.data, experiments=training.set,
    method = "L-BFGS-B", control=list(maxit=500, factr=1e9, trace=1, REPORT=1),
    lower = fit.lower[names(p)],  upper = fit.upper[names(p)])
  fit.par[names(fit$par)] <- fit$par
  cat(paste0('fit.par = c(', paste0(names(fit.par), ' = ', round(fit.par,4), collapse = ', '), ')'))
}

#############################################
#      Run Model
#############################################
for (exp.index in 1:length(experiments)) {
  exp.data = all.data[exp==experiments[exp.index]]
  somic.out = as.data.frame(somic(exp.data,
                              mic_vmax = fit.par['mic_vmax'],
                              mic_km = fit.par['mic_km'],
                              kdissolution = fit.par['kdissolution'],
                              kdepoly = fit.par['kdepoly'],
                              kdeath_and_exudates = fit.par['kdeath_and_exudates'],
                              kdesorb = fit.par['kdesorb'],
                              ksorb = fit.par['ksorb'],
                              kmicrobial_uptake = fit.par['kmicrobial_uptake'],
                              cue_0 = fit.par['cue_0'],
                              mcue = fit.par['mcue'],
                              mclay = fit.par['mclay'],
                              clay = exp.const$clay[exp.index], 
                              max_tsmd = exp.const$max_tsmd[exp.index],
                              use_atsmd = TRUE))
  suppressWarnings(
    all.data[exp==experiments[exp.index], (names(somic.out)) := somic.out]
  )
}

all.data[, facet_label := gsub('expt_','',exp)]
all.data[, facet_label := gsub('_',' ',facet_label)]
all.data[, facet_label := gsub('pendleton','pen',facet_label)]
all.data[, facet_label := gsub('sanborn','san',facet_label)]
all.data[exp %in% training.set, facet_label := paste(facet_label,'*')]

plot.smooth = FALSE
if (plot.smooth) {
  all.data[, soc.smth := regSmooth(time, soc, lambda = 3e-4), by = exp]
} else {
  all.data[, soc.smth := soc]
}

# insert some dummy panels to leave blank...
# ...as position holders to be replaced with regression plot
blank_panels = data.table(facet_label = c('bb22a', 'bb22b', 'bb5a', 'bb5b', 'hoos6a', 'hoos6b'))
all.data = rbindlist(list(all.data, blank_panels), fill = TRUE)

longterm.plots = ggplot(all.data, aes(x=time/12)) +
  geom_line(aes(y=soc.smth, color='Predicted'),size = 0.4, show.legend =FALSE) +
  geom_point(aes(y=measured.soc, shape = 'Observed'), data=all.data[!is.na(measured.soc)],
    color="black", size=2, show.legend=FALSE) +
  scale_color_manual(name = "", values = 'black', breaks = 'Predicted') +
  scale_shape_discrete(name = "", solid = FALSE, breaks = 'Observed') +
  scale_y_continuous(name = expression(SOC~'(Mg ha'^{-1}*')'), limits=c(0,90)) +
  scale_x_continuous(name = 'Time (years)') +
  facet_wrap(~facet_label, ncol=4) +
  theme_bw() +
  theme(panel.spacing = unit(0, "lines"),
    panel.border = element_rect(colour = 'grey50', fill=NA, size=0.3),
    strip.text = element_text(size=10),
    axis.text = element_text(size=8),
    axis.title =element_text(size=10),
    legend.position = c(0.95, 0.9), legend.justification = c(0.9, 0),
    legend.spacing = unit(0.1,'cm'), legend.margin = margin(-6,0,-6,0),
    legend.background = element_rect(color='transparent'),
    legend.box.background =  element_rect(color='transparent'))

modelled.soc   = all.data[!is.na(measured.soc) & !(exp %in% training.set), soc]
measured.soc   = all.data[!is.na(measured.soc) & !(exp %in% training.set), measured.soc]
mod.fit <- lm(modelled.soc~measured.soc)
summary(mod.fit)

somic.fit.data = data.table(
  modelled = modelled.soc,
  measured = measured.soc,
  model = 'SOMic')

lm_eqn = function(df){
  m = lm(modelled ~ measured, df)
  eq = substitute(
    atop(italic(y) == b * italic(x) + a,~~italic(r)^2~"="~r2),
    list(a  = format(coef(m)[1], digits = 2, nsmall=2),
         b  = format(coef(m)[2], digits = 2, nsmall=2),
         r2 = format(summary(m)$r.squared, digits = 2, nsmall=2)))
  as.character(as.expression(eq))
}
eq = somic.fit.data[, .(lm_eqn(data.frame(modelled,measured)))]

regress.plot = ggplot(somic.fit.data, aes(measured, modelled)) +
  geom_abline(slope = 1, color = 'grey80', size = 0.3) +
  geom_point(shape=21, stroke=0.2) +
  geom_smooth(method = "lm", se=FALSE, color="firebrick", linetype = 'longdash', formula = y ~ x, size = 0.8) +
  geom_text(data=eq, aes(40, 90, label=V1), parse = TRUE, inherit.aes=FALSE, size=3) +
  scale_x_continuous(name = expression(Observed~SOC~'(Mg ha'^{-1}*')'), limits = c(0,100)) +
  scale_y_continuous(name = expression(Predicted~SOC~'(Mg ha'^{-1}*')'), limits = c(0,100)) +
  coord_equal() +
  theme(axis.text = element_text(size=10),
    axis.title =element_text(size=10))

############### all plots with regression inset ################### 
# library(ggplotify)
library(grid)
g <- plot_to_gtable(longterm.plots)
# delete the empty panels tp provide room for regression plot
td = which(g$layout$l > 9 & g$layout$t <20)
g$grobs[td] <- NULL
g$layout <- g$layout[-td, ]

setEPS()
postscript(paste0(results_path,"LongTermPlots.eps"), width = 6, height = 6)
grid.newpage()
grid.draw(g)
vp = viewport(x=.53, y=.58, width=.42, height=.42, just = c('left','bottom'))
pushViewport(vp)
grid.draw(plot_to_gtable(regress.plot))
upViewport()
dev.off()


