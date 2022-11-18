



mainOutputFolder <- paste0('output/outputs_', as.character(Sys.time(), format = "%Y%m%d-%H%M%S"))

#
dir.create(mainOutputFolder, recursive = T)
file.copy('runOptimizations.R', file.path(mainOutputFolder, 'used-runOptimizations.R'))

nRepRandom <- 1024 # 1024
nRepBO <- 160 # 160

totalIterBO <- 20
propPointsBO <- 4
nCpuBO <- propPointsBO

seed_simSetup <- 2022
seed_randomSchemes <- 1877
seed_BOrandomMarkers <- 1809
seed_BOrealMarkers <- 1890
seed_BOestimatedMarkers <- 1822

outFolder_simSetup <- file.path(mainOutputFolder, 'simSetup')
outFolder_randomSchemes <- file.path(mainOutputFolder, 'randomSchemes')
outFolder_BOrandomMarkers <- file.path(mainOutputFolder, 'BOrandomMarkers')
outFolder_BOrealMarkers <- file.path(mainOutputFolder, 'BOrealMarkers')
outFolder_BOestimatedMarkers <- file.path(mainOutputFolder, 'BOestimatedMarkers')


meanTimeSimu <- 20
(expectedTimeRandom <- (meanTimeSimu * ceiling(nRepRandom / 64)) / 60)
(expectedTimeRandom <- (meanTimeSimu * ceiling((nRepBO * (totalIterBO + 1) + 1) / (64/nCpuBO))) /60)


# create simulation setup
rmarkdown::render('createSimulationSetup.Rmd',
                  params = list(
                    nGen = 5,
                    he = 0.7,
                    plotBudjetPerGen = 600,
                    outFolder = outFolder_simSetup,
                    seed = seed_simSetup),
                  output_dir = mainOutputFolder)
setupFile <- list.files(outFolder_simSetup, full.names = TRUE)


# TEST repetability: run random schemes
tryCatch({
  rmarkdown::render('runRandomBreedingSchemes.Rmd',
                    params = list(
                      nRepetition = 10,
                      simSetupFile = setupFile,
                      outFolder = file.path(mainOutputFolder, 'TEST-randomSchemes'),
                      seed = 1234),
                    output_file = file.path(mainOutputFolder,'test-randomSchemes'))
}, error = function(err) {
  message("Error detected for: runRandomBreedingSchemes.Rmd. Here's the original error message:")
  message(err)
  return(1)
})




# # run random schemes
# tryCatch({
#   rmarkdown::render('runRandomBreedingSchemes.Rmd',
#                     params = list(
#                       nRepetition = nRepRandom,
#                       simSetupFile = setupFile,
#                       outFolder = outFolder_randomSchemes,
#                       seed = seed_randomSchemes),
#                     output_dir = mainOutputFolder)
# }, error = function(err) {
#   message("Error detected for: runRandomBreedingSchemes.Rmd. Here's the original error message:")
#   message(err)
#   return(1)
# })





# run BO random markers effects
tryCatch({
  rmarkdown::render('runBO-randomMarkerEffects.Rmd',
                    params = list(
                      nRepetition = nRepBO,
                      totalIter =  totalIterBO,
                      proposedPoints = propPointsBO,
                      nCpuBO = nCpuBO,
                      simSetupFile = setupFile,
                      BOdirectory = file.path(outFolder_BOrandomMarkers, 'boResults'),
                      outFolder = outFolder_BOrandomMarkers,
                      seed = seed_BOrandomMarkers),
                    output_dir = mainOutputFolder)
}, error = function(err) {
  message("Error detected for: runBO-randomMarkerEffects.Rmd. Here's the original error message:")
  message(err)
  return(1)
})

# run BO estimated markers effects
tryCatch({
  rmarkdown::render('runBO-estimatedMarkerEffects.Rmd',
                    params = list(
                      nRepetition = nRepBO,
                      totalIter =  totalIterBO,
                      proposedPoints = propPointsBO,
                      nCpuBO = nCpuBO,
                      simSetupFile = setupFile,
                      BOdirectory = file.path(outFolder_BOestimatedMarkers, 'boResults'),
                      outFolder = outFolder_BOestimatedMarkers,
                      seed = seed_BOestimatedMarkers),
                    output_dir = mainOutputFolder)
}, error = function(err) {
  message("Error detected for: runBO-estimatedMarkerEffects.Rmd. Here's the original error message:")
  message(err)
  return(1)
})

# # run BO real markers effects
# tryCatch({
#   rmarkdown::render('runBO-RealMarkerEffects.Rmd',
#                     params = list(
#                       nRepetition = nRepBO,
#                       totalIter =  totalIterBO,
#                       proposedPoints = propPointsBO,
#                       nCpuBO = nCpuBO,
#                       simSetupFile = setupFile,
#                       BOdirectory = file.path(outFolder_BOrealMarkers, 'boResults'),
#                       outFolder = outFolder_BOrealMarkers,
#                       seed = seed_BOrealMarkers),
#                     output_dir = mainOutputFolder)
# }, error = function(err) {
#   message("Error detected for: runBO-RealMarkerEffects.Rmd. Here's the original error message:")
#   message(err)
#   return(1)
# })
