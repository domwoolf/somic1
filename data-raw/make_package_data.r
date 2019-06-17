time_step = 'month'
source("SOMIC packages.r")
source("SOMIC Functions.r")
exp.const = fread("experiment_constants.csv")
exp.const[, max_tsmd := -(20+1.3*clay - 0.01*clay^2) * depth/23]
raw.soc.data = get_data(".")
all.data = initialise.somic.data(copy(raw.soc.data))
usethis::use_data(exp.const, all.data, fit.par, overwrite = T)
