# This file demonstrates use of the SOMic model 
# using data from long-term experimental stations
library(data.table)
library(ggplot2)
data("exp.const")
data("all.data")
data("fit.par")
setDT(exp.const)
setDT(all.data)

experiments = unique(exp.const$exp)
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

ggplot(all.data, aes(x=time/12)) +
  geom_line(aes(y=soc, color='Predicted'),size = 0.4, show.legend =FALSE) +
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
