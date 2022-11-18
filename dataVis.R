
library(ggplot2)
library(dplyr)

dir_main <- 'output/1stTest'
dir_checkDtaCpu <- list.files(dir_main, pattern = 'TEST', full.names = TRUE)

dir_dta <- c(
  dir_randomSchemes = file.path(dir_main, 'randomSchemes'),
  dir_BOestim = file.path(dir_main, 'BOestimatedMarkers'),
  dir_BOrandom = file.path(dir_main, 'BOrandomMarkers'),
  dir_BOreal = file.path(dir_main, 'BOrealMarkers')
)



# Check repeatability between both servers ----
dta_checkCpu <- lapply(dir_checkDtaCpu, function(dir){
  dtaFile <- list.files(dir,
                        pattern = "aggregate",
                        full.names = TRUE)
  if (length(dtaFile) == 1) {
    dta <- read.csv(dtaFile)
    cols <- c("i", "iHomo", "bRep", "phenoFreq", "seed", "BV_mean",
              "method", "nGen", "he", "plotBudjetPerGen")
    return(dta[, cols])
  } else {
    warning("file not found")
    return(NULL)
  }
})
all(sapply(dta_checkCpu, function(dta){identical(dta, dta_checkCpu[[1]])}))




# load data from differents methods ----
dta_list <- lapply(dir_dta, function(dir){
  dtaFile <- list.files(dir,
                        pattern = "aggregate",
                        full.names = TRUE)
  if (length(dtaFile) == 1) {
    return(read.csv(dtaFile))
  } else {
    warning("file not found")
    return(NULL)
  }
})

dta <- do.call(bind_rows, dta_list)
dta <- dta[order(dta$BV_mean, decreasing = T),]
dta <- group_by(dta, method)

dta_summary <- summarise(dta, meanBv = mean(BV_mean), varBv = var(BV_mean))


# t.test(x = filter(dta, method == 'BO-randomEffects')$BV_mean,
#        y = filter(dta, method == 'BO-estimatedEffects')$BV_mean)
#
# lm_mod <- lm('BV_mean ~ method', data = dta)
# summary(lm_mod)
# anova(lm_mod)
#
# aov_mod <- aov(BV_mean ~ 1 + method, data = dta)
# summary(aov_mod)
# TukeyHSD(aov_mod)


# Data visualisation: ----
set.seed(1234)
ggplot(dta, aes(x = method, y = BV_mean,  col = method)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter()

ggplot(dta, aes(x = method, y = BV_mean,  col = method)) +
  geom_boxplot()

# Basic density
ggplot(dta, aes(x=BV_mean, col = method)) +
  geom_density()

