
library(ggplot2)
library(dplyr)
library(purrr)
library(tidyr)
library(dplyr)
library(FactoMineR)
library(factoextra)

dir_main <- 'output/scenario1'
# dir_main <- 'output/scenario2'
dir_main <- c('output/scenario1', 'output/scenario2')


all_dta <- lapply(dir_main, function(dir_main){
  dir_checkDtaCpu <- list.files(dir_main, pattern = 'TEST', full.names = TRUE)

  dir_dta <- c(
    dir_randomSchemes = file.path(dir_main, 'randomSchemes'),
    dir_BOestim = file.path(dir_main, 'BOestimatedMarkers'),
    # dir_BOestim = file.path(dir_main, 'BOestimatedMarkers-2'),
    dir_BOrandom = file.path(dir_main, 'BOrandomMarkers'),
    dir_BOreal = file.path(dir_main, 'BOrealMarkers'),
    dir_BOsuccesive = file.path(dir_main, 'BOsuccessiveOpt')
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


  # BV of the initial population
  simSetupDir <- paste0(dir_main, '/simSetup')
  simSetupFile <- list.files(simSetupDir, '.rds$', full.names = TRUE)
  simSetup <- readRDS(simSetupFile)
  gvInitPop <- simSetup$fixedParams$trait$gv(simSetup$fixedParams$initPop)
  summary(gvInitPop)
  (meanGvInitPop <- mean(gvInitPop))
  gvInitPop <- data.frame(
    mean = meanGvInitPop,
    scenario = paste(
      simSetup$fixedParams$nGen,
      simSetup$setupInfo$he,
      simSetup$fixedParams$plotBudjetPerGen,
      sep = "/"
    )
  )


  # load data from differents methods ----
  dta_list <- lapply(dir_dta, function(dir){
    dtaFiles <- list.files(dir,
                           pattern = ".csv$",
                           full.names = TRUE)
    # aggregateFile <- dtaFiles[which(grepl("aggregate", dtaFiles))]
    dtaFiles <- dtaFiles[!grepl('aggregate', dtaFiles)]
    aggregateFile <- c()
    if (length(dtaFiles) == 0) {
      warning("No data file found")
      return(NULL)
    } else if (length(aggregateFile) == 1) {
      out <- read.csv(aggregateFile)
    } else {
      out <- do.call(bind_rows,lapply(dtaFiles, read.csv))
    }
    out$i <- as.character(out$i)
    out$bRep <- as.character(out$bRep)
    out$phenoFreq <- as.character(out$phenoFreq)
    out$file <- as.character(dtaFiles)
    return(out)
  })
  dta <- do.call(bind_rows, dta_list)

  list(dta = dta,
       meanGvInitPop = gvInitPop)
})
all_dta <- all_dta %>% transpose()
dta <- do.call(bind_rows, all_dta$dta)
dta <- dta[order(dta$BV_mean, decreasing = T),]
dta$scenario <- paste(dta$nGen, dta$he, dta$plotBudjetPerGen, sep = "/")
dta$iInit <- dta$iHomo

dta$method[dta$method == "BO-realEffects"] <- "Actual QTN"
dta$methodOrder[dta$method == "Actual QTN"] <- 5
dta$method[dta$method == "BO-estimatedEffects"] <- 'Estimated QTNs effect'
dta$methodOrder[dta$method == "Estimated QTNs effect"] <- 2
dta$method[dta$method == "BO-randomEffects"] <- 'Random QTNs effect'
dta$methodOrder[dta$method == "Random QTNs effect"] <- 1
dta$method[dta$method == "BO-successiveOpt"] <- 'Repeated optimisation'
dta$methodOrder[dta$method == "Repeated optimisation"] <- 3
dta$method[dta$method == "random"] <- 'Random Schemes'
dta$methodOrder[dta$method == "Random Schemes"] <- 4

dta$method <- reorder(dta$method, dta$methodOrder) # reorder factor levels

dta$scenario[dta$scenario == "10/0.3/200"] <- "H2 = 0.3; nGen = 10; B/gen = 200"
dta$scenario[dta$scenario == "5/0.7/600"] <- "H2 = 0.7; nGen = 5; B/gen = 600"

meanGvInitPop <- do.call(bind_rows, all_dta$meanGvInitPop)
meanGvInitPop$scenario[meanGvInitPop$scenario == "10/0.3/200"] <- "H2 = 0.3; nGen = 10; B/gen = 200"
meanGvInitPop$scenario[meanGvInitPop$scenario == "5/0.7/600"] <- "H2 = 0.7; nGen = 5; B/gen = 600"




# get simulation parameters
source('src/simulation.R')
source('src/simulationFromGen2.R')
source('src/simulationWithInitPheno.R')
simParams <- do.call(bind_rows, lapply(1:nrow(dta), function(l){
  x <- dta[l,]
  optParams <- lapply(x[,c('i', "iHomo", "bRep", "phenoFreq")], function(p) {
    if (is.character(p)) {
    as.numeric(unlist(strsplit(x = p, split = '/')))
    } else{p}
  })
  remainBudget <- as.numeric(x['plotBudjetPerGen']) * as.numeric(x['nGen'])

  if (x$method %in% c('Random Schemes',
                      'Random QTNs effect',
                      'Actual QTN')) {
    simulParams <- getSimulParams(i = optParams$i[1],
                                  iHomo = optParams$iHomo[1],
                                  bRep = optParams$bRep[1],
                                  phenoFreq = optParams$phenoFreq[1],
                                  budget = remainBudget,
                                  nGen = as.numeric(x['nGen']),
                                  nIndIni = 198,
                                  plotCost = 1,
                                  newIndCost = 1)
  }
  if (x$method == 'Estimated QTNs effect') {
    simulParams <- getSimulParams2(i = optParams$i[1],
                                   iHomo = optParams$iHomo[1],
                                   bRep = optParams$bRep[1],
                                   phenoFreq = optParams$phenoFreq[1],
                                   budget = remainBudget,
                                   nIniPheno = 2*198,
                                   nGen = as.numeric(x['nGen']),
                                   nIndIni = 198,
                                   plotCost = 1,
                                   newIndCost = 1)
  }
  if (x$method == 'Repeated optimisation') {
    simulParams <- getSimulParams(i = optParams$i[1],
                                  iHomo = optParams$iHomo[1],
                                  bRep = optParams$bRep[1],
                                  phenoFreq = optParams$phenoFreq[1],
                                  budget = remainBudget,
                                  nGen = as.numeric(x['nGen']),
                                  nIndIni = 198,
                                  plotCost = 1,
                                  newIndCost = 1)
    remainBudget <- remainBudget - simulParams$nPheno[1] - simulParams$nNew[1]
    simulParams <- data.frame(as.list(sapply(simulParams,"[[",1))) # get first value of each list element


    for (gen in seq_along(optParams$i)[-1]) {
      if (remainBudget < 0) {
        browser()
      }
      p <- getSimulParams3(i = optParams$i[gen],
                           bRep =  optParams$bRep[gen],
                           phenoFreq = optParams$phenoFreq[gen],
                           budget = remainBudget,
                           nGen = as.numeric(x['nGen']),
                           nIndIni = simulParams[gen-1, 'nNew'],
                           plotCost = 1,
                           newIndCost = 1,
                           currentSelectionCycle = gen)
      remainBudget <- remainBudget - p$nPheno[1] - p$nNew[1]
      simulParams[gen,] <- sapply(p,"[[",1)
    }
    simulParams <- as.list(simulParams)
    simulParams$eff.budget <- sum(simulParams$nPheno) + sum(simulParams$nNew)
    simulParams$eff.i <- simulParams$nSelected/c(198, simulParams$nNew[-length(simulParams$nNew)])
  }

  simulParams <- as.data.frame(simulParams)
  simulParams$gen <- 1:nrow(simulParams)

  simulParams <- pivot_wider(simulParams,
                           names_from = gen,
                           values_from = c(nPheno, nSelected, nNew, eff.i))
  simulParams
}))

dta <- bind_cols(dta, simParams)






# View(filter(dta, method == "BO-successiveOpt"))
# View(filter(dta, method == "BO-randomEffects"))


dta <- group_by(dta, method, scenario)
(dta_summary <- summarise(dta, n = n(),  meanBv = mean(BV_mean), varBv = var(BV_mean), t = mean(simulationTime)))


dtaSummary <- filter(dta, method == "Random Schemes", nGen == 10)
dtaSummary$i <- as.numeric(dtaSummary$i)
dtaSummary$bRep <- as.numeric(dtaSummary$bRep)
dtaSummary$phenoFreq <- as.numeric(dtaSummary$phenoFreq)
summary(dtaSummary)

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
theme_set(theme_bw())
theme_update(text = element_text(size = 10))
outDir <- 'figures'

ggsave <- function(file, ...) {
  file <- paste0(file, '.jpg') # .svg
  ggplot2::ggsave(filename = file, width = 20, height = 15, units = "cm", dpi = 300)
}

set.seed(1234)
dta <- dta[order(dta$BV_mean),]
(
  p <- ggplot(dta, aes(x = method, y = BV_mean, col = method)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha=I(0.4)) +
    geom_hline(data = meanGvInitPop, aes(yintercept = mean, group = scenario)) + #, yintercept = mean, linetype="dashed", color = scenario) +
    facet_wrap(~scenario, scale = "fixed") +
    labs(y = "Mean breeding value of the final population") +
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
)


file <- file.path(outDir, 'BO2-resultsBoxplot')
suppressWarnings(dir.create(dirname(file), recursive = TRUE))
ggsave(file, p)
p


ggplot(dta, aes(x = method, y = BV_mean, col = method)) +
  geom_boxplot() +
  geom_hline(data = meanGvInitPop, aes(yintercept = mean, group = scenario)) +
  facet_wrap(~scenario, scale = "fixed")

# Basic density
ggplot(dta, aes(x = BV_mean, col = method)) +
  geom_density() +
  facet_wrap(~scenario, scale = "fixed") # +
# geom_vline(data = meanGvInitPop, aes(xintercept = mean, group = scenario))












(
  p <- ggplot(dta, aes(x = method, y = nPheno_1, col = method)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(alpha=I(0.4)) +
    geom_hline(data = meanGvInitPop, aes(yintercept = mean, group = scenario)) + #, yintercept = mean, linetype="dashed", color = scenario) +
    facet_wrap(~scenario, scale = "fixed") +
    labs(y = "Mean breeding value of the final population") +
    theme(axis.text.x = element_text(angle = -45, vjust = 1, hjust=0))
)


dtaPlot <- filter(dta, method %in% c('Actual QTN', 'Random QTNs effect', 'Repeated optimisation'))
pheno1 <- data.frame(
  v = rep(396, 2),
  scenario = unique(dta$scenario)
)

(p <- ggplot(dtaPlot, aes(x = nPheno_1, col = method)) +
    geom_density() +
    geom_vline(data = pheno1, aes(xintercept = v, group = scenario)) + #, yintercept = mean, linetype="dashed", color = scenario) +
    xlim(0, 1500) +
    facet_wrap(~scenario, scale = "fixed")
)

file <- file.path(outDir, paste0('BO2-nPheno1-density'))
suppressWarnings(dir.create(dirname(file), recursive = TRUE))
ggsave(file, p)
p



## PCA on the opt params of each method

library(FactoMineR)

# PCA
subdta <- filter(dta, scenario == "H2 = 0.3; nGen = 10; B/gen = 200")
subdta <- filter(dta, scenario == "H2 = 0.7; nGen = 5; B/gen = 600")
# subdta <- filter(subdta, method != c("Repeated optimisation"))
subdta <- filter(subdta, method != c("Random Schemes"))
subdta <- ungroup(subdta)
lapply(dtaList_scenario, function(subdta){

  scenario <- unique(subdta$scenario)
  columns <- c(
    which(colnames(subdta) == "method"),
    which(grepl('nPheno_', colnames(subdta))),
    which(grepl('nSelected_', colnames(subdta))),
    which(grepl('nNew_', colnames(subdta)))
  )
  pcaDta <- subdta[, columns]
  pcaDta <- pcaDta[,which(unlist(lapply(pcaDta, function(x)!all(is.na(x))))),with=F]
  # pcaDta$i <- as.numeric(pcaDta$i)
  # pcaDta$iInit <- as.numeric(pcaDta$iInit)
  # pcaDta$bRep <- as.numeric(pcaDta$bRep)
  # pcaDta$phenoFreq <- as.numeric(pcaDta$phenoFreq)



  pcaRes <- PCA(pcaDta,
                scale.unit = TRUE,
                # scale.unit = FALSE,
                quali.sup = 1,
                ncp = 4,
                graph = T)


  # add PCs to the data
  subdta <- mutate(subdta,
                   pca1 = pcaRes$ind$coord[, 1],
                   pca2 = pcaRes$ind$coord[, 2])

  scenarioGravityCenter <- as.data.frame(pcaRes$quali.sup$coord)
  scenarioGravityCenter$method <- row.names(scenarioGravityCenter)
  scenarioGravityCenter$method <- factor(scenarioGravityCenter$method, levels = levels(subdta$method))
  # scenarioGravityCenter$method <- reorder(scenarioGravityCenter$method, levels(subdta$method))
 str(scenarioGravityCenter$method)


 str(subdta$method)
  class(scenarioGravityCenter$method)

  plotDta <- subdta
  # plotDta <- filter(subdta, scenario %in% c('a: he_0.3_nGen_5_b_200\n',
  #                                        'h: he_0.7_nGen_10_b_600\n'))
  plotGravity <- scenarioGravityCenter
  # plotGravity <- filter(scenarioGravityCenter,
  #                       scenario %in% c('a: he_0.3_nGen_5_b_200\n',
  #                                       'h: he_0.7_nGen_10_b_600\n'))
  p <- ggplot(plotDta, aes(x = pca1,
                           y = pca2,
                           col = method))
  p <- p + geom_vline(xintercept = 0)
  p <- p + geom_hline(yintercept = 0)
  p <- p + geom_density_2d(aes())
  p <- p + geom_point(alpha = 4/10) # add all points
  # add variables arrows
  p <- p + geom_point(data = plotGravity,
                      size = 4,
                      shape=23, color="black",
                      aes(x = Dim.1,
                          y = Dim.2,
                          fill = method))
  p <- p + labs(title = scenario,
                fill = 'Center of gravity')
  # p <- p + annotate("segment",
  #                   x = 0,
  #                   xend = pcaRes$var$coord[,1],
  #                   y = 0,
  #                   yend = pcaRes$var$coord[,2],
  #                   size = 0.5,
  #                   arrow = arrow(length = unit(0.2, "cm")))
  # p <- p + annotate("text",
  #                   x = pcaRes$var$coord[,1],
  #                   y = pcaRes$var$coord[,2],
  #                   label = row.names(pcaRes$var$coord),
  #                   hjust = 'outward',
  #                   vjust = 'outward')
  axisLim <- c(min(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[1],
                   ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1])-0.1,
               max(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[2],
                   ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[2]+0.1)
  )
  # p <- p + xlim(axisLim) + ylim(axisLim)
  # suppressMessages({
  #   p <- p + coord_fixed(ratio = 1,
  #                        xlim = axisLim,
  #                        ylim = axisLim)
  # })
  # suppressMessages({
  #   p <- p + coord_fixed(ratio = 1,
  #                        xlim = c(ggplot_build(p)$layout$panel_scales_x[[1]]$range$range[1],
  #                                 5),
  #                        ylim = c(ggplot_build(p)$layout$panel_scales_y[[1]]$range$range[1],
  #                                 5))
  # })
  suppressMessages({
    p <- p + coord_fixed(ratio = 1,
                         xlim = c(axisLim[1]+2,# +2
                                  6+2),
                         ylim = c(axisLim[1],
                                  6))
  })
  p

  file <- file.path(outDir, paste0('BO2-resultsPCA_nGen-', unique(subdta$nGen)))
  suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  ggsave(file, p)
  p




  p <- (fviz_pca_var(pcaRes, repel = T)
        + labs(
          title = scenario)
        # '; optimization repetition: ',
        # unique(dtaList[[m]]$id)))
        + theme(plot.background = element_rect(fill = "white", color = "white")))
  file <- file.path(outDir, paste0('BO2-resultsPCA_Variables_nGen-', unique(subdta$nGen)))
  suppressWarnings(dir.create(dirname(file), recursive = TRUE))
  ggsave(file, p)
  p


})









source('src/simulationFromGen2.R')
source('src/simulation.R')
dtaSucc <- filter(dta, method == 'BO-successiveOpt')
simParams <- apply(dtaSucc, MARGIN = 1, function(x) {
  # browser()

  optParams <- lapply(x[c('i', "bRep", "phenoFreq")], function(p) {
    as.numeric(unlist(strsplit(x = p, split = '/')))
  })
  optParams$iHomo <- as.numeric(x['iHomo'])

  remainBudget <- as.numeric(x['plotBudjetPerGen']) * as.numeric(x['nGen'])
  gen1 <- getSimulParams(i = optParams$i[1],
                         iHomo = optParams$iHomo[1],
                         bRep = optParams$bRep[1],
                         phenoFreq = optParams$phenoFreq[1],
                         budget = remainBudget,
                         nGen = as.numeric(x['nGen']),
                         nIndIni = 198,
                         plotCost = 1,
                         newIndCost = 1)

  remainBudget <- remainBudget - gen1$nPheno - gen1$nNew

  simParams <- data.frame(as.list(sapply(gen1,"[[",1))) # get first value of each list element


    for (gen in seq_along(optParams$i)[-1]) {
      p <- getSimulParams3(i = optParams$i[gen],
                           bRep =  optParams$bRep[gen],
                           phenoFreq = optParams$phenoFreq[gen],
                           budget = remainBudget,
                           nGen = as.numeric(x['nGen']),
                           nIndIni = simParams[gen-1, 'nNew'],
                           plotCost = 1,
                           newIndCost = 1,
                           currentSelectionCycle = gen)
      remainBudget <- remainBudget - p$nPheno - p$nNew

      simParams[gen,] <- sapply(p,"[[",1)
    }
  simParams[,'gen'] <- 1:nrow(simParams)
  simParams[,'method'] <- x['method']
  simParams[,'scenario'] <- x['scenario']
  simParams[,'mainSeed'] <- x['mainSeed']
  simParams
})


dtaBORand <- filter(dta, method == 'BO-randomEffects')

simParams <- c(simParams, apply(dtaBORand, MARGIN = 1, function(x) {

  simParams <- getSimulParams(i = as.numeric(x['i']),
                              iHomo = as.numeric(x['iHomo']),
                              bRep =  as.numeric(x['bRep']),
                              phenoFreq = as.numeric(x['phenoFreq']),
                              budget = as.numeric(x['plotBudjetPerGen']) * as.numeric(x['nGen']),
                              nGen = as.numeric(x['nGen']),
                              nIndIni = 198,
                              plotCost = 1,
                              newIndCost = 1)
  simParams <- as.data.frame(simParams)
  simParams[,'gen'] <- 1:nrow(simParams)
  simParams[,'method'] <- x['method']
  simParams[,'scenario'] <- x['scenario']
  simParams[,'mainSeed'] <- x['mainSeed']
  simParams
}))

simParams <- do.call(rbind, simParams)
simParams$bRep <- simParams$nNew/(simParams$nPheno + simParams$nNew)


seeds <- simParams[c(min(which(simParams$method == 'BO-randomEffects')), min(which(simParams$method == 'BO-successiveOpt'))), 'mainSeed']

simParams2 <- simParams[simParams$mainSeed %in% seeds,]

simParams_s1 <- filter(simParams, scenario == "5/0.7/600")
ggplot(simParams_s1, aes(x=gen, y=nNew, col=method, group = mainSeed)) +
  geom_point(size=2) +
  geom_line()

simParams_s1$gen <- as.factor(simParams_s1$gen)
ggplot(simParams_s1, aes(x = gen, y=nNew, col = method)) +
  geom_boxplot()




simParams_s2 <- filter(simParams, scenario == "10/0.3/200")
ggplot(simParams_s2, aes(x=gen, y=nNew, col=method, group = mainSeed)) +
  geom_point(size=2) +
  geom_line()

simParams_s2$gen <- as.factor(simParams_s2$gen)
ggplot(simParams_s2, aes(x = gen, y=nNew, col = method)) +
  geom_boxplot() +
  coord_cartesian(ylim = c(0, 500))



s <- as.numeric(simParams_s1[simParams_s1$nNew == max(simParams_s1$nNew), 'mainSeed'][1])
simParams_s1[as.numeric(simParams_s1$mainSeed) == s,]










dtaEstS2 <- filter(dta, method == "Estimated QTNs effect")
dtaEstS2 <- filter(dtaEstS2, scenario == "H2 = 0.3; nGen = 10; B/gen = 200")

boSeeds <- sort(dtaEstS2$BOseed)
boSeeds2 <- sort(as.numeric(regmatches(list.files('output/scenario2/BOestimatedMarkers/boResults/'),
           regexpr("\\d+",
                   list.files('output/scenario2/BOestimatedMarkers/boResults/')))))
boSeeds == boSeeds2

worstScheme <- filter(dta, method == "Estimated QTNs effect")
worstScheme <- worstScheme[worstScheme$BV_mean %in% sort(worstScheme$BV_mean)[1:3],]
worstScheme[, c('iHomo', 'i', 'phenoFreq', 'bRep')]
source('src/simulationWithInitPheno.R')
s <- 1


worstScheme[s, c('iHomo', 'i', 'phenoFreq', 'bRep', "BV_mean")]
getSimulParams2(i = as.numeric(worstScheme$i[s]),
                iHomo = as.numeric(worstScheme$iHomo[s]),
                bRep = as.numeric(worstScheme$bRep[s]),
                phenoFreq = as.numeric(worstScheme$phenoFreq[s]),
                budget = as.numeric(worstScheme$plotBudjetPerGen[s]) * as.numeric(worstScheme$nGen[s]),
                nGen = as.numeric(worstScheme$nGen[s]),
                nIniPheno = 198*2,
                plotCost = 1,
                nIndIni = 198,
                newIndCost = 1)
worstScheme$BOseed
worstScheme$method
optFile <- paste0('output/scenario2/BOestimatedMarkers/boResults/boRes_noName__', worstScheme$BOseed[s], '.rds')
worstOpt <- readRDS(optFile)
           paste0('output/scenario2/BOestimatedMarkers/boResults/boRes_noName__', bestScheme$BOseed[s], '.rds')
optFile2 <- paste0('output/scenario2/BOestimatedMarkers/boResults/boRes_noName__100873.rds')
otherOpt <- readRDS(optFile2)

View(worstOpt$results[,
                      c("mu", "nQTN",
                        "nCpus", "nCpusFunEval",
                        "totalIter",
                        "i",
                        "iHomo", "bRep", "phenoFreq", "BV_mean", "dob", "eol", "error.message",
                        "exec.time", "ei", "error.model", "train.time", "prop.type",
                        "propose.time", "se", "mean", "seed")
])

library(xtable)

plotWorstDta <- worstOpt$results[,c("dob", "BV_mean", "i",
                                    "iHomo", "bRep", "phenoFreq")]
names(plotWorstDta)[1] <- "Iteration"
xtable(plotWorstDta)

p <- (ggplot(plotWorstDta, aes(x = Iteration, y = BV_mean))
      + geom_point()
)
file <- file.path(outDir, paste0('BO2-worstOpt'))
suppressWarnings(dir.create(dirname(file), recursive = TRUE))
ggsave(file, p)

#hist
# plot_ly(x = worstOpt$simSetup$fixedParams$trait$qtnEff,
#         type = "histogram",
#         histnorm = "",# "" | "percent" | "probability"
#         name = "dta",
#         # nbinsx = 100, # no more than nbinsx
#         hovertemplate = "x: (%{x})\ny: %{y}",
#         marker = list(color = "green",
#                       line = list(color = "black",
#                                   width = 1))
# ) %>%
#   layout(
#     title = "x",
#     xaxis = list(title = "x"),
#     yaxis = list(title = "Count", # "percent" | "probability"
#                  zeroline = TRUE)
#   )

hist(worstOpt$simSetup$fixedParams$trait$qtnEff)
hist(otherOpt$simSetup$fixedParams$trait$qtnEff)


length(worstOpt$simSetup$fixedParams$trait$qtnEff)
worstOpt$simSetup$fixedParams$trait$gv(worstOpt$simSetup$fixedParams$initPop)
pheno_W <- worstOpt$simSetup$fixedParams$initPhenoData

actualSetup <- readRDS('output/scenario1/simSetup/simSetup_f4ccf7249703731dfa9822ae4dc17b74.rds')

#plot
plot_ly(type = "scatter",
        mode = "lines+markers",
        x = actualSetup$fixedParams$trait$gv(actualSetup$fixedParams$initPop),
        y = worstOpt$simSetup$fixedParams$trait$gv(worstOpt$simSetup$fixedParams$initPop)
)

actGV <- actualSetup$fixedParams$trait$gv(actualSetup$fixedParams$initPop)
worstGV <- worstOpt$simSetup$fixedParams$trait$gv(worstOpt$simSetup$fixedParams$initPop)
bestGV <- bestOpt$simSetup$fixedParams$trait$gv(bestOpt$simSetup$fixedParams$initPop)

plot(x=actGV,
     y=worstGV)
cor(x=actGV,
     y=worstGV)
plot(x=actGV,
     y= bestGV)
cor(x=actGV,
     y=bestGV)
min()
max()


bestScheme <- dtaEstS2[dtaEstS2$BV_mean %in% sort(dtaEstS2$BV_mean, decreasing = T)[1:3],]
source('src/simulationWithInitPheno.R')
bestScheme[, c('iHomo', 'i', 'phenoFreq', 'bRep')]


s <- 2
getSimulParams2(i = as.numeric(bestScheme$i[s]),
                iHomo = as.numeric(bestScheme$iHomo[s]),
                bRep = as.numeric(bestScheme$bRep[s]),
                phenoFreq = as.numeric(bestScheme$phenoFreq[s]),
                budget = as.numeric(bestScheme$plotBudjetPerGen[s]) * as.numeric(bestScheme$nGen[s]),
                nGen = as.numeric(bestScheme$nGen[s]),
                nIniPheno = 198*2,
                plotCost = 1,
                nIndIni = 198,
                newIndCost = 1)
bestScheme$BOseed
bestScheme$method
optFile <- paste0('output/scenario2/BOestimatedMarkers/boResults/boRes_noName__', bestScheme$BOseed[s],'.rds')
bestOpt <- readRDS(optFile)
bestOpt$simSetup$fixedParams$trait$gv(bestOpt$simSetup$fixedParams$initPop)
pheno_B <- bestOpt$simSetup$fixedParams$initPhenoData

bestOpt$simSetup$fixedParams$phenotyper$ve
bestOpt$simSetup$fixedParams$phenotyper$mu
worstOpt$simSetup$fixedParams$phenotyper$ve
worstOpt$simSetup$fixedParams$phenotyper$mu

pheno <- cbind(pheno_B, pheno_W)

#plot
plot_ly(type = "scatter",
        mode = "markers",
        x = pheno_B$trait1,
        y = pheno_W$trait1
)
c
