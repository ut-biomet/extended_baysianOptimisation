# Author: Julien Diot juliendiot@ut-biomet.org
# 2020 The University of Tokyo
#
# Description:
# Definition of function setting up the simulation

#' Setup the simulation for the optimization
#'
#' @param dataFile genetic data of the initial population. ("vcf" file)
#' @param lchrCm length of the chromosome in centi-morgan
#' @param nQTN total number of QTN
#' @param nGen Number of generation
#' @param plotCost cost for phenotyping one plot
#' @param newIndCost cost for creating one new individual
#' @param selectMateInds selection and mating function
#' @param createModel function for creating prediction model
#' @param breedSimOpt Breeding simulation function parametrized with optimized
#' and fixed parameters
#' @param aggrFun function aggregating the calculated BV
#' @param aggrFunName
#' @param totalBudget total budget of the breeding campaign
#' @param plotBudgetPerGen budget express as number of plotCost per generation
#' the total budget will be calculated by:
#' `totalBudget = nGen * plotCost * plotBudgetPerGen`
#' @param nSNP total number of SNP
#' @param genoChipSize number of SNP for genotyping
#' @param genoChipNind number of individuals to subset for selecting the
#' genotyped SNP
#' @param genoChipMinMAF threshold on the minor allele frequency in the subset
#' of individuals. Only SNP with MAF >= genoChipMinMAF will be selected (randomly).
#' @param mu phenotypic mean (see: breedSimulatR::phenotyper)
#' @param he heritability (see: breedSimulatR::phenotyper)
#' @param ve environmental variance (see: breedSimulatR::phenotyper)
#' @param specName specie's name
#' @param popName population name
#' @param traitName trait name
#' @param outputFolder folder where to save the setup informations. Set it to
#' NULL to not save the results. Default "simSetups".
#' @param replace If they already exists, should the setup files be replaced or
#' saved with a new name. (default: TRUE)
#' @param setupName Name of this setup. It will be recored in the setup, and
#' will be used to create the file's name of the saved file. Default: "noName"
#' @param seed [optional] random seed
#' @param verbose [optional] Display information
#' @param genotypes [optional] genetic data of the initial population (bed matrix)
#' @param qtnEff [optional] Effects of the QTNs
#'
#' @return list of setup information: specie,snps,initPop = initPop,trait and
#' phenotyper. (see: breedSimulatR's package documentation)
setupSimulation3 <- function(dataFile=NULL,
                             lchrCm,
                             nQTN,
                             nGen,
                             plotCost,
                             newIndCost,
                             selectMateInds,
                             createModel,
                             breedSimOpt,
                             aggrFun,
                             aggrFunName,
                             totalBudget = NULL,
                             plotBudjetPerGen = NULL,
                             nSNP = NULL,
                             genoChipSize = NULL,
                             genoChipNind = 20,
                             genoChipMinMAF = 0.1,
                             currentSelectionCycle,
                             currentModel,
                             genoDta,
                             phenoDta,
                             currentPop,
                             mu = 0,
                             he = NULL,
                             ve = NULL,
                             specName = "bayesOpt",
                             popName = "initialPop",
                             traitName = "trait1",
                             outputFolder = 'simSetups',
                             replace = TRUE,
                             setupName = "noName",
                             varSetupParams = NULL,
                             seed = NULL,
                             genotypes = NULL,
                             qtnEff = NULL,
                             verbose = TRUE) {

  # test totalBudget and plotBudjetPerGen
  if ((is.null(totalBudget) && is.null(plotBudjetPerGen))
    || (!is.null(totalBudget) && !is.null(plotBudjetPerGen))) {
    stop('`totalBudget` or (exlusif) `plotBudjetPerGen` should be specified')
  }


  # set random seed:
  if (!is.null(seed)) {
    prevSeed <- .Random.seed
    set.seed(seed)
  }


  # Create breedSimulatR's objects ----
  BSR_Obj <- createBreedSimObj3(dataFile = dataFile,
                                specName = specName,
                                lchrCm = lchrCm,
                                popName = popName,
                                nQTN = nQTN,
                                traitName = traitName,
                                mu = mu,
                                ve = ve,
                                he = he,
                                qtnEff = currentModel$qtnEff,
                                verbose = verbose)



  # Select SNP for genotyping ----
  if (is.null(genoChipSize)) {
    genoChipSNP <- colnames(BSR_Obj$initPop$genoMat)
  } else {
    genoChipPop <- BSR_Obj$initPop$clone()
    # select a subset of individuals from the population
    genoChipPop$remInds(indsNames = sample(names(genoChipPop$inds),
                                           genoChipPop$nInd - genoChipNind))

    # sample markers for which the MAF in the samples individuals is above
    # specific threshold
    maf <- genoChipPop$maf
    genoChipSNP <- names(maf[maf >= genoChipMinMAF])

    if (length(genoChipSNP) >= genoChipSize) {
      genoChipSNP <- sample(genoChipSNP, genoChipSize)
    } else {
      stop("Couldn't sample SNP for genotyping.",
           "Not enought markers after filtering on the MAF.")
    }
  }



  # Create fixed parameters list ----

  # calculate total budget
  if (is.null(totalBudget)) {
    totalBudget <- nGen * plotCost * plotBudjetPerGen
  }

  fp <- list() # list of fixed parameters (constraints)
  fp$nGen           = nGen
  fp$budget         = totalBudget
  fp$plotBudjetPerGen = totalBudget / nGen
  fp$plotCost       = plotCost
  fp$newIndCost     = newIndCost
  fp$initPop        = BSR_Obj$initPop      # initial population
  fp$trait          = BSR_Obj$trait        # trait of interest
  fp$phenotyper     = BSR_Obj$phenotyper   # phenotyper
  fp$genoChipSNP    = genoChipSNP
  fp$selectMateInds = selectMateInds
  fp$createModel    = createModel
  fp$aggrFun        = aggrFun
  fp$aggrFunName    = aggrFunName
  fp$currentSelectionCycle = currentSelectionCycle
  fp$currentModel = currentModel
  fp$genoDta = genoDta
  fp$phenoDta = phenoDta
  fp$currentPop = currentPop

  initPopId        = paste0("initPop", "_", digest(BSR_Obj$initPop))
  traitId          = paste0("trait", "_", digest(BSR_Obj$trait))
  phenotyperId     = paste0("phenoLab", "_", digest(BSR_Obj$phenotyper))
  createModelId    = paste0("createModel", "_", digest(createModel))
  selectMateIndsId = paste0("selectMateInds", "_", digest(selectMateInds))
  aggrFunNameId    = paste0(aggrFunName, "_", digest(aggrFun))

  # Create optimized parameters: ----
  minInds <- 3
  minBrep <- (minInds * fp$newIndCost) / fp$budget
  maxBrep <- 1

  optParams <- makeParamSet(
    makeNumericParam(id = "i",
                     lower = 0.01,
                     upper = 0.99,
                     default = 0.5),
    makeNumericParam(id = "bRep",
                     lower = minBrep,
                     upper = maxBrep,
                     default = 0.5),
    makeIntegerParam(id = "phenoFreq",
                     lower = 1,
                     upper = fp$nGen,
                     default = 1)
  )



  # Create Objective function ----
  objFun <- createObjFun3(fixedParams = fp,
                          breedSimOpt = breedSimOpt)


  # Finalization ----
  setup <- list(
    setupInfo = list(dataFile = dataFile,
                     dataFileId = digest::digest(file = dataFile),
                     nSNP = nSNP,
                     setupSeed = seed,
                     filteredDataId = NULL,
                     nQTN = nQTN,
                     mu = mu,
                     ve = ve,
                     he = he,
                     initPopId = initPopId,
                     traitId = traitId,
                     phenotyperId = phenotyperId,
                     createModelId = createModelId,
                     selectMateIndsId = selectMateIndsId,
                     aggrFunNameId = aggrFunNameId),
    setupName = setupName,
    fixedParams = fp,
    objFun = objFun,
    optParams = optParams,
    varOfInterest = varSetupParams
  )

  setupHash <- digest::digest(setup)
  if (verbose) {
    cat(paste("Simulation setup MD5 fingerprint:", setupHash, "\n"))
  }


  # save setup in RDS file ----
  if (!is.null(outputFolder)) {
    resultFile <- paste0(outputFolder,
                         "/simSetup_",
                         setupName)
    if (is.null(setupName)) {
      resultFile <- paste0(resultFile, setupHash)
    }
    resF <-  paste0(resultFile, '.rds')
    if (!replace) {
      i <- 1
      while (file.exists(resF)) {
        i <- i + 1
        resF <- paste0(resultFile, '_', i, '.rds')
      }
    }
    saveRDS(setup, file = resF)

    if (verbose) {
      fileHash <- digest::digest(resF, file = TRUE)
      cat(paste("Setup file MD5 fingerprint:", fileHash, "\n"))
      cat('File created:', resF, '\n')
    }
  }


  # reset the seed
  if (!is.null(seed)) {
    set.seed(prevSeed)
  }

  setup
}


createBreedSimObj3 <- function(dataFile,
                               specName,
                               lchrCm,
                               nChr,
                               popName,
                               nQTN = NULL,
                               traitName,
                               mu,
                               ve,
                               he,
                               qtnEff = NULL,
                               verbose) {


  ## specie definition ---------------------


  if (verbose) cat("\nCreate breedSimulatR's pop, snps and specie objects\n")
  breeSimObj <- breedSimulatR::readVCF(file = dataFile, specie = NULL, verbose = verbose)
  breeSimObj$pop

  ## trait definition ---------------------
  if (verbose) cat("\nCreate breedSimulatR's trait object\n")
  if (is.null(nQTN)) {
    nQTN <- breeSimObj$snps$nSNP()
  } else if ((nQTN <= 0) | (nQTN > breeSimObj$snps$nSNP())) {
    stop("nSNP should be a positive number lower or equal to nSNP or to the total number of SNP in the provided data.")
  }

  ## sample qtn
  qtnId <- sample(breeSimObj$snps$SNPcoord$SNPid, nQTN)


  ## qtn effects
  if (is.null(qtnEff)) {
    varQtnEffects <- 6 / (2 * sum(breeSimObj$pop$af[qtnId] * (1 - breeSimObj$pop$af[qtnId])))
    qtnEffects <- rexp(nQTN, sqrt(2/(varQtnEffects))) * sample(c(-1, 1), nQTN, replace = T)
    names(qtnEffects) <- qtnId
  } else {
    qtnEffects <- qtnEff
  }

  trait1 <- breedSimulatR::trait$new(name = traitName,
                                     qtn = qtnId,
                                     qtnEff = qtnEffects[qtnId])

  ## phenotyper definition ----------------
  if (verbose) cat("\nCreate breedSimulatR's phenotyper object\n")
  if (!is.null(ve )) {
    phenoLab <- breedSimulatR::phenotyper$new(name = "PhenoLab1",
                                              traits = list(trait1),
                                              plotCost = 1,
                                              mu = mu,
                                              ve = ve,
                                              he = NULL,
                                              pop = NULL)
  } else {
    phenoLab <- breedSimulatR::phenotyper$new(name = "PhenoLab1",
                                              traits = list(trait1),
                                              plotCost = 1,
                                              mu = mu,
                                              ve = ve,
                                              he = he,
                                              pop = breeSimObj$pop)
  }

  ## output ----
 list(initPop = breeSimObj$pop,
      trait = trait1,
      phenotyper = phenoLab)
}








createObjFun3 <- function(fixedParams,
                          breedSimOpt) {

  # Create a simple function ----
  # This function is doing the simulation with only one vector (x)
  # as parameter. This is necessary for using `mlrMBO`:
  fun <- function(x) {
    breedSimOpt(i = x[1],
                bRep = x[2],
                phenoFreq = x[3],
                seed = x[4],
                budget = fixedParams$budget,
                nGen = fixedParams$nGen,
                currentSelectionCycle = fixedParams$currentSelectionCycle,
                currentModel = fixedParams$currentModel,
                genoDta = fixedParams$genoDta,
                phenoDta = fixedParams$phenoDta,
                currentPop = fixedParams$currentPop,
                plotCost = fixedParams$plotCost,
                newIndCost = fixedParams$newIndCost,
                trait = fixedParams$trait,
                phenotyper = fixedParams$phenotyper,
                genoChipSNP = fixedParams$genoChipSNP,
                createModel = fixedParams$createModel,
                selectMateInds = fixedParams$selectMateInds,
                aggrFun = fixedParams$aggrFun,
                aggrFunName = fixedParams$aggrFunName,
                verbose = FALSE)
  }
  fun
}
