# Title: haploRILs_function.R
# Author: jmontero
# Date: 2024-10-11
# Description: Required functions for the execution of haploRILs 

#### Currently used functions ####

#### SNP SCORE RULES ####
# When a SNP matches the founders on
# a minor allele, it gets the score
minimum_af_allele_match = 1
# When a SNP matches the founders on
# a major allele, it gets the score
maximum_af_allele_match = 1
# When a SNP does not match the founders
# or there is missing data
no_match_or_missing_data = 0

#### BUILD VARIABLE-SIZE WINDOWS WITH MINIMUM BP AND SNP ####
# Step: 0 --> Stable window, no sliding
# Step: 4 --> slide jump = 4
getSlidingWindows = function(map = map, minimum_snps, minimum_distance = 0, step) {
  # Initialize coord vectors
  binStart_bp = c()
  binEnd_bp = c()
  binStart_snps = c()
  binEnd_snps = c()
  # Initialize index and iterate
  i = 1
  while (i < (nrow(map))) {
    # Add start coordinates
    binStart_bp = append(binStart_bp, map$pos[i])
    binStart_snps = append(binStart_snps, i)
    # Update the window size for this window, based on the minimum distance criterium
    variableWindowSize = minimum_snps - 1
    while ((map$pos[i+variableWindowSize] - map$pos[i]) < minimum_distance  && (i + variableWindowSize) <= nrow(map)) {
      variableWindowSize = variableWindowSize + 1
    }
    # Add end coordinates (considering end of data frame)
    if (((i+variableWindowSize) <= nrow(map)) & (((i+variableWindowSize) <= (nrow(map)-minimum_snps)) |
                                                 ((map$pos[nrow(map)]-map$pos[i+variableWindowSize]) >= minimum_distance))) {
      binEnd_bp = append(binEnd_bp, map$pos[i+variableWindowSize])
      binEnd_snps = append(binEnd_snps, i+variableWindowSize)
    } else {
      binEnd_bp = append(binEnd_bp, map$pos[nrow(map)])
      binEnd_snps = append(binEnd_snps, nrow(map))
      break
    }
    # Update the index based on the step
    #i = ifelse(step == 0, i + variableWindowSize + 1, i + variableWindowSize/step + 1)
    if (step == 0) {
      i = i + variableWindowSize + 1
    } else {
      i = i + step
    }
  }
  # Make final table
  binCoords =
    data.frame(bin_MP = round((binStart_bp + binEnd_bp)/2, 0),
               binStart_bp = binStart_bp,
               binEnd_bp = binEnd_bp,
               binStart_snps = binStart_snps,
               binEnd_snps = binEnd_snps,
               binlength_bp = binEnd_bp - binStart_bp,
               binlength_snps = binEnd_snps - binStart_snps+1,
               density_snps_kbpsnps = (binEnd_snps - binStart_snps+1)/((binEnd_bp - binStart_bp)/1000))
  return(binCoords)
}

#### DETERMINE THE MINIMUM SNP NUMBER PER WINDOW TO REACH THE AIM RESOLUTION ####
# Determine window step (fixed window step)
defineMinimumNumberSNPs = function(map, aimResolution = 2e3, noLessThan_SNPs = 5) {
  # Chromosome length
  chromLength = map$pos[nrow(map)]
  # Number of loci
  snpNumber = nrow(map)
  # How many windows to achieve the target resolution
  numberWindows = chromLength/aimResolution
  # Average number of SNPs per window
  windowSize = round(snpNumber/numberWindows, 0)
  # Is window size under 5?
  condition = windowSize < noLessThan_SNPs
  # If window size is under 5, increase the resolution
  while ( condition ) {
    # Add 500 bp
    aimResolution = aimResolution + 0.5e3
    # Update condition
    windowSize = round(snpNumber/(chromLength/aimResolution), 0)
    condition = windowSize < noLessThan_SNPs
  }
  # Chromosome SNP density
  density = (chromLength/1000)/snpNumber
  # Minimum distance based on SNP density
  minimum_distance = 1000*density*windowSize
  # Determine window step (fixed window step)
  step = round(windowSize/4, 0)
  print(
    paste0("chrom length: ", chromLength, 
           " bp - number SNPs: ", nrow(map) ,
           " - resolution: ", aimResolution,
           " bp - number windows: ", numberWindows,
           " SNPs - window size: ", windowSize, 
           " SNPs - minimum distance: ", minimum_distance, 
           " bp - step: ", step, " SNPs")
  )
  return(windowSize)
}

#### GET STATS OF DATA FRAMES ####
getStats = function(df) {
  print(deparse(substitute(df)))
  report =
    df %>%
    summarise(no_windows = nrow(df),
              avr_SNPs = mean(binlength_snps),
              avr_bp = mean(binlength_bp))
  print(report)
}

#### BUILD SLIDING WINDOWS ####
# Maybe writing a function for building sliding
# windows is overkill, but it complies with the
# previous function getBinCoords() that required
# a bit more coding
getBinCoords = function(map, minimum_snps, minimum_distance) {
  binStart_bp = c()
  binEnd_bp = c()
  binStart_snps = c()
  binEnd_snps = c()
  n = 1
  while (n < (nrow(map))) {
    binStart_bp = append(binStart_bp, map$pos[n])
    binStart_snps = append(binStart_snps, n)
    e = minimum_snps - 1
    while ((map$pos[n+e] - map$pos[n]) < minimum_distance  && (n + e) <= nrow(map)) {
      e = e + 1
    }
    if (((n+e) <= nrow(map)) & (((n+e) <= (nrow(map)-minimum_snps)) |
                                ((map$pos[nrow(map)]-map$pos[n+e]) >= minimum_distance))) {
      binEnd_bp = append(binEnd_bp, map$pos[n+e])
      binEnd_snps = append(binEnd_snps, n+e)
    } else {
      binEnd_bp = append(binEnd_bp, map$pos[nrow(map)])
      binEnd_snps = append(binEnd_snps, nrow(map))
      break
    }
    n = n + e + 1
  }
  binCoords =
    data.frame(bin_MP = round((binStart_bp + binEnd_bp)/2, 0),
               binStart_bp = binStart_bp,
               binEnd_bp = binEnd_bp,
               binStart_snps = binStart_snps,
               binEnd_snps = binEnd_snps,
               binlength_bp = binEnd_bp - binStart_bp,
               binlength_snps = binEnd_snps - binStart_snps+1,
               density_snps_kbpsnps = (binEnd_snps - binStart_snps+1)/((binEnd_bp - binStart_bp)/1000))
  return(binCoords)
}


#### ASSIGN SCORES TO FOUNDER SNPS BASED ON SIMILARITY WITH THE DESCENDANT ####
# The set number 2 should reward founders that have less the
# less common genotypes matching with the RILs" genotypes
# and therefore make clearer distinction between the
# best-scoring founders
assignScoresRules2 = function(ped) {
  ped = as.matrix(ped)
  scores = sapply(1:ncol(ped), function(j) {
    score_col = ifelse(ped[-nrow(ped), j] != ped[nrow(ped), j],
                       no_match_or_missing_data, ## No match --> Remove 1
                       ifelse(ped[-nrow(ped), j] == ped[nrow(ped), j] & ped[-nrow(ped), j] == 2,
                              maximum_af_allele_match, ## Match on major allele --> Add 1
                              ifelse(ped[-nrow(ped), j] == ped[nrow(ped), j] & ped[-nrow(ped), j] == 1,
                                     minimum_af_allele_match, ## Match on major allele --> Add 1
                                     no_match_or_missing_data)) ## Both genotypes are missing data
    )
    return(score_col)
  })
  colnames(scores) = colnames(ped)
  return(scores)
}

#### SUM SCORES BY WINDOWS AND DETECT HIGHEST-SCORING FOUNDERS ####
# blockWindowBuilder takes a matrix of founder scores, with founder scores, and
# window parameters and calculates the sum score by windows. Then, it records
# the founders with the maximum score and extracts the rownames. Optionally, it
# applies inferMostLikelyFounder() to the sum score matrix
blockWindowBuilder = function(scores, binCoords, genotypes = scores, corWindowLength, export = "mostLikelyFounders", exportMatrix = FALSE) {
  #### Get the by-window score sum ####
  scoreSum = matrix(nrow = nrow(scores), ncol = nrow(binCoords))
  binNames = c()
  for (w in 1:nrow(binCoords)) {
    scoresWindow = scores[, binCoords$binStart_snps[w]:binCoords$binEnd_snps[w]]
    binNames = append(binNames, paste(colnames(scoresWindow)[c(1, ncol(scoresWindow))], collapse = "-"))
    scoreSum[, w] = as.numeric(round(rowSums(scoresWindow)/binCoords$binlength_snps[w], 4))
  }
  rownames(scoreSum) = rownames(scores)
  colnames(scoreSum) = binNames
  #### Find highest-score founders ####
  maxScore = matrix(nrow = 1
                    , ncol = ncol(scoreSum))
  for (j in 1:ncol(scoreSum)) {
    mostLikelyFounder = paste(names(which(scoreSum[, j] == max(scoreSum[, j]))), collapse = "/")
    maxScore[, j] = mostLikelyFounder
  }
  colnames(maxScore) = binNames
  rownames(maxScore) = "mostLikelyFounder"
  #### Get matrix of highest-scoring founder names ####
  haplotypePaths = matrix(0, nrow = nrow(scoreSum), ncol = ncol(scoreSum))
  for (j in 1:ncol(scoreSum)) {
    haplotypePaths[which(scoreSum[, j] == max(scoreSum[, j])), j] = 1
  }
  rownames(haplotypePaths) = rownames(scores)
  colnames(haplotypePaths) = binNames
  #### Calculate most-likely haplotype sequence using the function calculate_most_likely_path ####
  bestHaplotypePaths = as.matrix(inferMostLikelyPathAlternative6(scoreSum, k = corWindowLength, slidingWindows = binCoords, genotypes = genotypes, export = "mostLikelyFounders"))
  #rownames(bestHaplotypePaths) = binNames
  colnames(bestHaplotypePaths) = "mostLikelyFounder"
  exportList = list(haplotypePaths = t(haplotypePaths),
                    maxScoreFounders = t(maxScore),
                    bestHaplotypePaths = bestHaplotypePaths)
  if (exportMatrix == FALSE) {
    return(bestHaplotypePaths)
  } else {
    return(scoreSum)
  }
  
}

#### FILTER HBLOCKS USING THE K-CONTEXT METHOD AND EXTRACT A TABLE WITH THE MOST-LIKELY FOUNDER PER HBLOCK ####
# There are different versions but haploRILs currently uses inferMostLikelyPathAlternative6
# Given a sequence of founder haplotypes, clean it up to have the most-likely
# haplotype sequence that assigns snps to single founders

#### inferMostLikelyPathAlternative6 ####
# Retained COs must be switches between non-recurrent haplotypes: none of the
# recombination interval haplotypes must match any of the context region ones
inferMostLikelyPathAlternative6 = function(df, k, slidingWindows, genotypes, export = "mostLikelyFounders") {
  # Tools for counting blocks
  blocks = c(1)
  block = 1
  # Tools for extracting putative CO intervals
  fineIntervals = c()
  fineBlocks = c()
  for (i in 2:ncol(df)) {
    #### Iterate over the founder sequence, stitch them up together into hblocks or detect putative crossovers ####
    # Compare the highest-score founders between the last
    # and the current position. If they are the same, keep
    # the same block, but if they differ, evaluate if the
    # recombination is true using the context
    foundersLastRow = which(df[, i-1] %in% max(df[, i-1]))
    foundersThisRow = which(df[, i] %in% max(df[, i]))
    commonFoundersRows = intersect(foundersLastRow, foundersThisRow)
    if (length(commonFoundersRows) > 0) {
      # If two consecutive SNPs have common founders --> no CO, same block
      blocks = append(blocks, block)
    } else {
      #### Validate recombination events using the K-context method ####
      # Confirm recombination events using the context
      # haplotype sequences
      # The context will be split into four regions: (1)
      # upwards, -1-k away from current haplotype, (2)
      # last haplotype, i-1, (3) current haplotype, i, (4)
      # downwards, +k away from current haplotype.
      # Comparing the upwards and the downwards context
      # regions will provide support or contradict the
      # putative crossover between the contrasting current
      # and next haplotypes
      upperLimit = ifelse(i-1-k >= 1, i-1-k, 1)
      bottomLimit = ifelse(i+k <= ncol(df), i+k, ncol(df))
      nextSNP = ifelse(i+1 <= ncol(df), i+1, ncol(df))
      upwardKwindow = as.matrix(df[, upperLimit:(i-2)])
      downwardKwindow = as.matrix(df[, nextSNP:bottomLimit])
      foundersUpwardKwindowRows = unlist(apply(upwardKwindow, 2, function(x) { which(x == max(x)) }))
      foundersDownwardKwindowRows = unlist(apply(downwardKwindow, 2, function(x) { which(x == max(x)) }))
      #             Putative crossover --> X
      ###############################################################################
      #  upward k window  #    last row    #   this row   #   downwards k window    #
      ###############################################################################
      noHaplotypeSwitch = any(foundersUpwardKwindowRows %in% foundersDownwardKwindowRows # Compare downwards and downwards
                              | foundersUpwardKwindowRows %in% foundersThisRow # Compare upwards with this SNP
                              | foundersDownwardKwindowRows %in% foundersLastRow) # Compare downwards with last SNP
      noHaplotypeSwitch = ifelse(k == 0, FALSE, noHaplotypeSwitch)
      # If any common founder is present both in the upwards
      # and downwards context regions, do not count as
      # block
      if (noHaplotypeSwitch) {
        blocks = append(blocks, block)
      } else {
        # If no common founder, count new block
        block = block + 1
        blocks = append(blocks, block)
        #### Extract crossover coordinates ####
        # Extract the genotypes within the start and end of
        # the putative CO interval detected
        candidateRegion = as.data.frame(
          t(genotypes[c(foundersLastRow[1], foundersThisRow[1]),
                      slidingWindows$binStart_snps[i-1]:slidingWindows$binEnd_snps[i]]))
        lastColNames = colnames(candidateRegion)
        colnames(candidateRegion) = c("upwardFounder", "downwardFounder")
        # Keep SNPs that are different between the two
        # recombining founders and then select the interval
        # between the first and the second
        candidateRegion$het = ifelse(candidateRegion$upwardFounder ==
                                       candidateRegion$downwardFounder,
                                     FALSE,
                                     TRUE)
        candidateRegion = candidateRegion %>% dplyr::filter(het)
        candidateRegion =
          candidateRegion %>%
          mutate(snp = rownames(candidateRegion),
                 hap = rleid(upwardFounder, downwardFounder))
        candidateRegion =
          suppressWarnings(suppressMessages(
            candidateRegion %>%
              group_by(hap) %>%
              slice(c(1, n())) %>%
              slice(c(2, 3)) %>%
              select(upwardFounder, snp, downwardFounder)
          ))
        fineInterval =
          paste0(lastColNames[1], ":", candidateRegion$snp[1], "-",
                 candidateRegion$snp[2], ":", lastColNames[2])
        fineIntervals = append(fineIntervals, fineInterval)
        fineBlocks = append(fineBlocks, block)
      }
    }
  }
  fineDF = data.frame(block = fineBlocks,
                      interval = fineIntervals)
  #### Create HBK output table, add interval coordinates ####
  # Create a table with the haploblock coordinates and most
  # likely founders
  blocksDF = data.frame(index = 1:length(blocks), block = blocks)
  blocksDF =
    suppressWarnings(suppressMessages(
      blocksDF %>%
        mutate(interval = colnames(df)[as.numeric(rownames(blocksDF))]) %>%
        group_by(block) %>%
        summarise(startIndex = first(index), endIndex = last(index),
                  startInterval = first(interval), endInterval = last(interval))
    ))
  blocksDF$startBP = as.numeric(str_match(blocksDF$startInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,3])
  blocksDF$startBP[1] = as.numeric(str_match(blocksDF$startInterval[1], "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,2])
  blocksDF$endBP = as.numeric(str_match(blocksDF$endInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,2])
  blocksDF$endBP[nrow(blocksDF)] = as.numeric(str_match(blocksDF$endInterval[nrow(blocksDF)], "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,3])
  blocksDF$startIntervalLength = as.numeric(str_match(blocksDF$startInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,3]) - as.numeric(str_match(blocksDF$startInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,2])
  blocksDF$endIntervalLength = as.numeric(str_match(blocksDF$endInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,3]) - as.numeric(str_match(blocksDF$endInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,2])
  #### Add additional cols about scores and most likely founders ####
  # Generate a new haplotype path that assigns 1 to the most
  # likely founders within each block
  mostLikelyFoundersBlocks = c()
  mostLikelyFoundersBlocks2 = c()
  newDF = df[]
  blockSumScoreMax1 = c()
  blockSumScoreMax2 = c()
  blockSumScoreRest = c()
  for (i in 1:nrow(blocksDF)) {
    rowSumsBlock = rowSums(as.matrix(df[, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)]))
    bestFounderRows = which(rowSumsBlock %in% max(rowSumsBlock))
    bestFounderRows2 = which(rowSumsBlock %in% unique(sort(rowSumsBlock, decreasing = TRUE))[2])
    mostLikelyFounders = paste(rownames(df)[bestFounderRows], collapse = "/")
    mostLikelyFoundersBlocks = append(mostLikelyFoundersBlocks, mostLikelyFounders)
    mostLikelyFounders2 = paste(rownames(df)[bestFounderRows2], collapse = "/")
    mostLikelyFoundersBlocks2 = append(mostLikelyFoundersBlocks2, mostLikelyFounders2)
    newDF[, blocksDF$startIndex[i]:blocksDF$endIndex[i]] = 0
    newDF[bestFounderRows, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)] = 1
    blockSumScoreMax1 = append(blockSumScoreMax1, mean(rowSumsBlock[bestFounderRows])/(blocksDF$endIndex[i] - blocksDF$startIndex[i] + 1))
    blockSumScoreMax2 = append(blockSumScoreMax2, mean(rowSumsBlock[bestFounderRows2])/(blocksDF$endIndex[i] - blocksDF$startIndex[i] + 1))
    blockSumScoreRest = append(blockSumScoreRest, mean(rowSumsBlock[-c(bestFounderRows, bestFounderRows2)])/(blocksDF$endIndex[i] - blocksDF$startIndex[i] + 1))
  }
  blocksDF$mostLikelyFounders = mostLikelyFoundersBlocks
  blocksDF$mostLikelyFounders2 = mostLikelyFoundersBlocks2
  blocksDF$blockSumScore1 = blockSumScoreMax1
  blocksDF$blockSumScore2 = blockSumScoreMax2
  blocksDF$blockSumScoreRest = blockSumScoreRest
  mostLikelyFoundersExport = c()
  for (i in 1:ncol(newDF)) {
    mostLikelyFoundersExport = append(mostLikelyFoundersExport,
                                      paste(rownames(newDF)[which(as.matrix(newDF[, i]) %in% 1)], collapse = "/"))
  }
  #### Export output, default = 'blocks' for HBK table ####
  exportList = list(as.data.frame(mostLikelyFoundersExport), fineDF, blocksDF)
  if (export == "mostLikelyFounders") {
    return (mostLikelyFoundersExport)
  } else if (export == "blocks") {
    return (blocksDF)
  } else if (export == "fineIntervals") {
    return (fineDF)
  } else if (export == "list") {
    return (exportList)
  }
}

#### IMPUTE REPETITIVE AND SOLATED DISSENTING BLOCKS ####
# Filter the blocks by merging repeating blocks and removing
# isolated dissenting blocks
filterBlocks = function(blocks) {
  newBlocks = data.frame()
  for (ID in unique(blocks$id)){
    # Merge consecutive haploblocks with same most likely
    # founder
    blocksID = subset(blocks, id == ID)
    blocksID = blocksID %>%
      select(mostLikelyFounders, startIndex, endIndex)
    blocksID = blocksID %>% mutate(groups = data.table::rleid(mostLikelyFounders))
    first_rows = blocksID %>% group_by(groups) %>% slice(1) %>% select(mostLikelyFounders, groups, startIndex)
    last_rows = blocksID %>% group_by(groups) %>% slice(n())  %>% select(endIndex, groups)
    blocksID = left_join(first_rows, last_rows, by = "groups")
    # Impute isolated dissenting haploblocks
    blocksID$lastLetter = dplyr::lag(blocksID$mostLikelyFounders)
    blocksID$nextLetter = lead(blocksID$mostLikelyFounders)
    blocksID =
      blocksID %>%
      mutate(newFounder = ifelse(lastLetter == nextLetter, lastLetter, mostLikelyFounders))
    blocksID[1, "newFounder"] = blocksID[1, "mostLikelyFounders"]
    blocksID[nrow(blocksID), "newFounder"] = blocksID[nrow(blocksID), "mostLikelyFounders"]
    blocksID = blocksID %>% transmute(mostLikelyFounders = newFounder, groups, startIndex, endIndex)
    blocksID$groups = rleid(blocksID$mostLikelyFounders)
    first_rows = blocksID %>% group_by(groups) %>% slice(1) %>% select(mostLikelyFounders, groups, startIndex)
    last_rows = blocksID %>% group_by(groups) %>% slice(n())  %>% select(endIndex, groups)
    blocksID = left_join(first_rows, last_rows, by = "groups")
    blocksID$id = ID
    blocksID = blocksID %>% select(id, everything())
    newBlocks = rbind(newBlocks, as.data.frame(blocksID))
  }
  blocks = as.data.frame(blocks) %>% mutate(startIndex = as.numeric(blocks$startIndex), endIndex = as.numeric(blocks$endIndex), step = as.numeric(blocks$step), K = as.numeric(blocks$K), nSnp = as.numeric(blocks$nSnp))
  newBlocks = as.data.frame(newBlocks) %>% mutate(id = as.character(id))
  newBlocks = as.data.frame(newBlocks) %>% rename(mostLikelyFoundersFiltered = mostLikelyFounders)
  blocks = left_join(as.data.frame(blocks),
                     as.data.frame(newBlocks %>% select(id, mostLikelyFoundersFiltered, startIndex, groups)),
                     by = c("id", "startIndex"))
  blocks = left_join(as.data.frame(blocks),
                     as.data.frame(newBlocks %>% select(id, mostLikelyFoundersFiltered, endIndex, groups)),
                     by = c("id", "endIndex"))
  blocks$blocksFiltered = as.numeric(dplyr::coalesce(blocks$groups.x, blocks$groups.y))
  blocks$mostLikelyFoundersFiltered = dplyr::coalesce(blocks$mostLikelyFoundersFiltered.x, blocks$mostLikelyFoundersFiltered.y)
  blocks = blocks %>% select(!ends_with(".x")) %>% select(!ends_with(".y"))
  blocks = blocks %>% tidyr::fill(mostLikelyFoundersFiltered) %>% tidyr::fill(blocksFiltered)
  blocks = blocks %>% select(id, chr, nSnp, step, K, blocksFiltered, block, startIndex, endIndex, startBP, endBP, mostLikelyFoundersFiltered, mostLikelyFounders, mostLikelyFounders2, blockSumScore1, blockSumScore2, blockSumScoreRest, startInterval, endInterval, startIntervalLength, endIntervalLength)
  intervalLengthDF =
    blocks %>%
    #group_by(id, nSnp, K, blocksFiltered) %>%
    group_by(id, nSnp, step, K, blocksFiltered, .groups = "drop") %>%
    slice(1) %>%
    summarise(startIntervalLength, endIntervalLength)
  intervalLengthDF$IntervalLength = (intervalLengthDF$startIntervalLength + intervalLengthDF$endIntervalLength)/2
  intervalLengthDF = intervalLengthDF %>% select(-startIntervalLength, -endIntervalLength)
  blocks = left_join(as.data.frame(blocks), as.data.frame(intervalLengthDF), by =c("id", "nSnp", "step", "K", "blocksFiltered")) %>% select(-startIntervalLength, -endIntervalLength)
  return(blocks)
}

#### Other functions ####
#### BUILD VARIABLE-SIZED SLIDING WINDOWS WITH MINIMUM BP AND SNP ####
getBinCoords = function(map, minimum_snps, minimum_distance) {
  binStart_bp = c()
  binEnd_bp = c()
  binStart_snps = c()
  binEnd_snps = c()
  n = 1
  while (n < (nrow(map))) {
    binStart_bp = append(binStart_bp, map$pos[n])
    binStart_snps = append(binStart_snps, n)
    e = minimum_snps - 1
    while ((map$pos[n+e] - map$pos[n]) < minimum_distance  && (n + e) <= nrow(map)) {
      e = e + 1
    }
    if (((n+e) <= nrow(map)) & (((n+e) <= (nrow(map)-minimum_snps)) |
                                ((map$pos[nrow(map)]-map$pos[n+e]) >= minimum_distance))) {
      binEnd_bp = append(binEnd_bp, map$pos[n+e])
      binEnd_snps = append(binEnd_snps, n+e)
    } else {
      binEnd_bp = append(binEnd_bp, map$pos[nrow(map)])
      binEnd_snps = append(binEnd_snps, nrow(map))
      break
    }
    n = n + e + 1
  }
  binCoords =
    data.frame(bin_MP = round((binStart_bp + binEnd_bp)/2, 0),
               binStart_bp = binStart_bp,
               binEnd_bp = binEnd_bp,
               binStart_snps = binStart_snps,
               binEnd_snps = binEnd_snps,
               binlength_bp = binEnd_bp - binStart_bp,
               binlength_snps = binEnd_snps - binStart_snps+1,
               density_snps_kbpsnps = (binEnd_snps - binStart_snps+1)/((binEnd_bp - binStart_bp)/1000))
  return(binCoords)
}


#### FILTER HBLOCKS USING THE K-CONTEXT METHOD AND EXTRACT A TABLE WITH THE MOST-LIKELY FOUNDER PER HBLOCK ####
# There are different versions but haploRILs currently uses inferMostLikelyPathAlternative6
# Given a sequence of founder haplotypes, clean it up to have the most-likely
# haplotype sequence that assigns snps to single founders

#### inferMostLikelyPath ####
# Basic algorithm, including forward and backward
# iterations, and basic context-based CO imputation by
# zeroing-out no-COs
inferMostLikelyPath = function(haplotypePaths, k, correctionMethod = FALSE) {
  df = haplotypePaths[]
  # First iteration: Forward. Keep only repetitive
  # haplotypes
  for (i in 2:ncol(df)) {
    bestFoundersThisRow = which(df[, i] == 1)
    bestFoundersLastRow = which(df[, i-1] == 1)
    commonFounders = intersect(bestFoundersThisRow, bestFoundersLastRow)
    if (length(commonFounders) > 0) {
      df[, i] = 0
      df[commonFounders, i] = 1
    }
  }
  # Second iteration: Backward. Keep repetitive haplotypes
  for (i in (ncol(df) - 1):1) {
    bestFoundersThisRow = which(df[, i] == 1)
    bestFoundersLastRow = which(df[, i+1] == 1)
    commonFounders = intersect(bestFoundersThisRow, bestFoundersLastRow)
    if (length(commonFounders) > 0) {
      df[, i] = 0
      df[commonFounders, i] = 1
    }
  }
  if (correctionMethod != FALSE) {
    # Third iteration: Confirm recombination events using the
    # context haplotype sequences
    for (i in (k+1):(ncol(df)-k-1)) {
      bestFoundersThisRow = which(df[, i] == 1)
      bestFoundersNextRow = which(df[, i+1] == 1)
      commonFounders = intersect(bestFoundersThisRow, bestFoundersNextRow)
      if (length(commonFounders) == 0 && sum(df[, i]) != 0) {
        upwardKwindow = df[, c((i-k):(i-1))]
        candidateWindow = df[, i:(i+1)]
        downwardKwindow = df[, c((i+1+1):(i+1+k))]
        upwardFounderRow = which(rowSums(upwardKwindow) %in% max(rowSums(upwardKwindow)))
        downwardFounderRow = which(rowSums(downwardKwindow) %in% max(rowSums(downwardKwindow)))
        identicalContextFoundersRow = intersect(upwardFounderRow, downwardFounderRow)
        if (length(identicalContextFoundersRow) != 0) {
          if (correctionMethod == "imputation") { # If imputation is selected, simply add 0 to these columns so that no founder can be annotated
            df[, i:(i+1)] = 0
          }
        }
      }
    }
  }
  # Fourth iteration: Extract the row name of the best
  # founder
  mostLikelyFounders= c()
  for (i in 1:ncol(df)) {
    mostLikelyFounders = append(mostLikelyFounders,
                                paste(rownames(df)[which(df[, i] == 1)], collapse = "/"))
  }
  return(mostLikelyFounders)
}

#### inferMostLikelyPathAlternative1 ####
# Modifies imputation technique. Now by copying over the
# last haplotype into the imputed gap
inferMostLikelyPathAlternative1 = function(haplotypePaths, k, correctionMethod = TRUE) {
  df = haplotypePaths[]
  # First iteration: Forward. Keep only repetitive
  # haplotypes
  for (i in 2:ncol(df)) {
    bestFoundersThisRow = which(df[, i] == 1)
    bestFoundersLastRow = which(df[, i-1] == 1)
    commonFounders = intersect(bestFoundersThisRow, bestFoundersLastRow)
    if (length(commonFounders) > 0) {
      df[, i] = 0
      df[commonFounders, i] = 1
    }
  }
  # Second iteration: Backward. Keep repetitive haplotypes
  for (i in (ncol(df) - 1):1) {
    bestFoundersThisRow = which(df[, i] == 1)
    bestFoundersLastRow = which(df[, i+1] == 1)
    commonFounders = intersect(bestFoundersThisRow, bestFoundersLastRow)
    if (length(commonFounders) > 0) {
      df[, i] = 0
      df[commonFounders, i] = 1
    }
  }
  if (correctionMethod != FALSE) {
    # Third iteration: Confirm recombination events using the
    # context haplotype sequences
    for (i in (k+1):(ncol(df)-k-1)) {
      bestFoundersThisRow = which(df[, i] == 1)
      bestFoundersNextRow = which(df[, i+1] == 1)
      commonFounders = intersect(bestFoundersThisRow, bestFoundersNextRow)
      if (length(commonFounders) == 0 && sum(df[, i]) != 0) {
        # The context will be split into four regions: (1)
        # upwards, -k away from current haplotype, (2)
        # current haplotype, i, (3) next haplotype, i+1, (4)
        # downwards, +1+k away from current haplotype.
        # Comparing the upwards and the downwards context
        # regions will provide support or contradict the
        # putative crossover between the contrasting current
        # and next haplotypes
        upwardKwindowAll = haplotypePaths[, c((i-k):(i-1))]
        downwardKwindowAll = haplotypePaths[, c((i+1+1):(i+1+k))]
        upwardKwindowAssigned = df[, c((i-k):(i-1))]
        downwardKwindowAssigned = df[, c((i+1+1):(i+1+k))]
        bestFounderUpwardKwindowRow = which(rowSums(upwardKwindowAssigned) %in% max(rowSums(upwardKwindowAssigned)))
        bestFounderDownwardKwindowRow = which(rowSums(downwardKwindowAssigned) %in% max(rowSums(downwardKwindowAssigned)))
        assignedFoundersUpwardKwindowRow = which(rowSums(upwardKwindowAll) %in% max(rowSums(upwardKwindowAll)))
        assignedFoundersDownwardKwindowRow = which(rowSums(downwardKwindowAll) %in% max(rowSums(downwardKwindowAll)))
        conditionNoHaplotypeSwitch = any(bestFounderUpwardKwindowRow %in% assignedFoundersDownwardKwindowRow) |
          any(bestFounderDownwardKwindowRow %in% assignedFoundersUpwardKwindowRow)
        # If the best founder is not the same between the
        # upward region and the current snp, impute the
        # current snp with the founder from the upward
        # region. This will discard false haploblocks that
        # survive imputation because their location within a
        # region with contrasting haplotypes
        if (!any(bestFounderUpwardKwindowRow %in% bestFoundersThisRow)) {
          df[, i] == 0
          df[bestFounderUpwardKwindowRow, i] == 1
        }
        # Do same for downward region
        if (!any(bestFounderDownwardKwindowRow %in% bestFoundersNextRow)) {
          df[, i+1] == 0
          df[bestFounderDownwardKwindowRow, i] == 1
        }
        # If the best founder upwards is present among the
        # assigned founders downwards or viceversa, impute
        # the haplotypes in the candidate crossover
        # interval, as the context haplotype info suggests
        # no haplotype switch. This is probably a
        # genotyping/alignment error
        if (conditionNoHaplotypeSwitch) {
          if (correctionMethod != FALSE) {
            df[, i:(i+1)] = 0
            df[bestFounderUpwardKwindowRow, i:(i+1)] = 1
          }
        }
      }
    }
  }
  # Fourth iteration: Extract the row name of the best
  # founder
  mostLikelyFounders= c()
  for (i in 1:ncol(df)) {
    mostLikelyFounders = append(mostLikelyFounders,
                                paste(rownames(df)[which(df[, i] == 1)], collapse = "/"))
  }
  return(mostLikelyFounders)
}

#### inferMostLikelyPathAlternative2 ####
# Compared to alt1, all founders are compared between the
# context regions. Also, considers k=0 as no imputation,
# which caused unexpected results in last algs
inferMostLikelyPathAlternative2 = function(df, k) {
  blocks = c(1)
  block = 1
  for (i in 2:ncol(df)) {
    # Compare the highest-score founders between the last
    # and the current position. If they are the same, keep
    # the same block, but if they differ, evaluate if the
    # recombination is true using the context
    foundersLastRow = which(df[, i-1] %in% 1)
    foundersThisRow = which(df[, i] %in% 1)
    commonFoundersRows = intersect(foundersLastRow, foundersThisRow)
    if (length(commonFoundersRows) > 0) {
      blocks = append(blocks, block)
    } else {
      # Confirm recombination events using the context
      # haplotype sequences
      
      # The context will be split into four regions: (1)
      # upwards, -1-k away from current haplotype, (2)
      # last haplotype, i-1, (3) current haplotype, i, (4)
      # downwards, +k away from current haplotype.
      # Comparing the upwards and the downwards context
      # regions will provide support or contradict the
      # putative crossover between the contrasting current
      # and next haplotypes
      upperLimit = ifelse(i-1-k >= 1, i-1-k, 1)
      bottomLimit = ifelse(i+k <= ncol(df), i+k, ncol(df))
      nextSNP = ifelse(i+1 <= ncol(df), i+1, ncol(df))
      upwardKwindow = as.matrix(df[upperLimit:(i-2)])
      downwardKwindow = as.matrix(df[nextSNP:bottomLimit])
      foundersUpwardKwindowRows = which(rowSums(upwardKwindow) %in% max(rowSums(upwardKwindow)))
      foundersDownwardKwindowRows = which(rowSums(downwardKwindow) %in% max(rowSums(downwardKwindow)))
      noHaplotypeSwitch = any(foundersUpwardKwindowRows %in% foundersDownwardKwindowRows)
      noHaplotypeSwitch = ifelse(k == 0, FALSE, noHaplotypeSwitch)
      
      # If any founder is present both in the upwards and
      # downwards context regions, do not count as block
      if (noHaplotypeSwitch) {
        blocks = append(blocks, block)
      } else {
        block = block + 1
        blocks = append(blocks, block)
      }
    }
  }
  # Create a table with the haploblock coordinates and most
  # likely founders
  blocksDF = data.frame(block = blocks)
  blocksDF = blocksDF %>%
    mutate(index = as.numeric(rownames(blocksDF))) %>%
    mutate(interval = colnames(df)[as.numeric(rownames(blocksDF))]) %>%
    group_by(block) %>%
    summarise(startIndex = first(index), endIndex = last(index),
              startInterval = first(interval), endInterval = last(interval))
  # Generate a new haplotype path that assigns 1 to the most
  # likely founders within each block
  mostLikelyFoundersBlocks = c()
  newDF = df[]
  for (i in 1:nrow(blocksDF)) {
    bestFounderRows =
      which(rowSums(as.matrix(df[, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)]))
            %in% max(rowSums(as.matrix(df[, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)]))))
    mostLikelyFounderBlocks = paste(rownames(df)[bestFounderRows], collapse = "/")
    mostLikelyFoundersBlocks = append(mostLikelyFoundersBlocks, mostLikelyFounderBlocks)
    newDF[, blocksDF$startIndex[i]:blocksDF$endIndex[i]] = 0
    newDF[bestFounderRows, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)] = 1
  }
  blocksDF$mostLikelyFounders = mostLikelyFoundersBlocks
  mostLikelyFoundersExport = c()
  for (i in 1:ncol(newDF)) {
    mostLikelyFoundersExport = append(mostLikelyFoundersExport,
                                      paste(rownames(df)[which(newDF[, i] %in% 1)],
                                            collapse = "/"))
  }
  return(mostLikelyFoundersExport)
}

#### inferMostLikelyPathAlternative3 ####
# Input is scoreSum matrix instead of haplotypePath.
# scoreSum contains the founder sum-scores per interval, so
# it might favour discrimination between founders by
# providing more accurate estimations. Notice that the
# haplotypePath gives 1 to the max-scoring founders and 0 to
# the others, so info is lost
inferMostLikelyPathAlternative3 = function(df, k, scoreThreshold) {
  blocks = c(1)
  block = 1
  for (i in 2:ncol(df)) {
    # Compare the highest-score founders between the last
    # and the current position. If they are the same, keep
    # the same block, but if they differ, evaluate if the
    # recombination is true using the context
    foundersLastRow = which(df[, i-1] >= scoreThreshold)
    foundersThisRow = which(df[, i] >= scoreThreshold)
    commonFoundersRows = intersect(foundersLastRow, foundersThisRow)
    if (length(commonFoundersRows) > 0 || length(foundersLastRow) != 0 || length(foundersThisRow) != 0) {
      blocks = append(blocks, block)
    } else {
      # Confirm recombination events using the context
      # haplotype sequences
      
      # The context will be split into four regions: (1)
      # upwards, -1-k away from current haplotype, (2)
      # last haplotype, i-1, (3) current haplotype, i, (4)
      # downwards, +k away from current haplotype.
      # Comparing the upwards and the downwards context
      # regions will provide support or contradict the
      # putative crossover between the contrasting current
      # and next haplotypes
      upperLimit = ifelse(i-1-k >= 1, i-1-k, 1)
      bottomLimit = ifelse(i+k <= ncol(df), i+k, ncol(df))
      nextSNP = ifelse(i+1 <= ncol(df), i+1, ncol(df))
      upwardKwindow = as.matrix(df[upperLimit:(i-2)])
      downwardKwindow = as.matrix(df[nextSNP:bottomLimit])
      foundersUpwardKwindowRows = which(rowSums(upwardKwindow) %in% max(rowSums(upwardKwindow)))
      foundersDownwardKwindowRows = which(rowSums(downwardKwindow) %in% max(rowSums(downwardKwindow)))
      noHaplotypeSwitch = any(foundersUpwardKwindowRows %in% foundersDownwardKwindowRows)
      noHaplotypeSwitch = ifelse(k == 0, FALSE, noHaplotypeSwitch)
      
      # If any founder is present both in the upwards and
      # downwards context regions, do not count as block
      if (noHaplotypeSwitch) {
        blocks = append(blocks, block)
      } else {
        block = block + 1
        blocks = append(blocks, block)
      }
    }
  }
  # Create a table with the haploblock coordinates and most
  # likely founders
  blocksDF = data.frame(block = blocks)
  blocksDF = blocksDF %>%
    mutate(index = as.numeric(rownames(blocksDF))) %>%
    mutate(interval = colnames(df)[as.numeric(rownames(blocksDF))]) %>%
    group_by(block) %>%
    summarise(startIndex = first(index), endIndex = last(index),
              startInterval = first(interval), endInterval = last(interval))
  # Generate a new haplotype path that assigns 1 to the most
  # likely founders within each block
  mostLikelyFoundersBlocks = c()
  newDF = df[]
  for (i in 1:nrow(blocksDF)) {
    bestFounderRows =
      which(rowSums(as.matrix(df[, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)]))
            %in% max(rowSums(as.matrix(df[, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)]))))
    mostLikelyFounderBlocks = paste(rownames(df)[bestFounderRows], collapse = "/")
    mostLikelyFoundersBlocks = append(mostLikelyFoundersBlocks, mostLikelyFounderBlocks)
    newDF[, blocksDF$startIndex[i]:blocksDF$endIndex[i]] = 0
    newDF[bestFounderRows, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)] = 1
  }
  blocksDF$mostLikelyFounders = mostLikelyFoundersBlocks
  mostLikelyFoundersExport = c()
  for (i in 1:ncol(newDF)) {
    mostLikelyFoundersExport = append(mostLikelyFoundersExport,
                                      paste(rownames(df)[which(newDF[, i] %in% 1)],
                                            collapse = "/"))
  }
  return(mostLikelyFoundersExport)
}

#### inferMostLikelyPathAlternative4 ####
# Adds the functionality to extract base pairs that are
# potentially linked to crossovers. This is after confirming
# crossovers, the genotypes of the bordering founders are
# compared and heterozygous snps extracted (match in one
# founder not in the other). Then, consecutive different
# heterozygous markers are saved
inferMostLikelyPathAlternative4 = function(df, k, slidingWindows, genotypes, export = "mostLikelyFounders") {
  # Tools for counting blocks
  blocks = c(1)
  block = 1
  # Tools for extracting putative CO intervals
  fineIntervals = c()
  fineBlocks = c()
  for (i in 2:ncol(df)) {
    # Compare the highest-score founders between the last
    # and the current position. If they are the same, keep
    # the same block, but if they differ, evaluate if the
    # recombination is true using the context
    foundersLastRow = which(df[, i-1] %in% max(df[, i-1]))
    foundersThisRow = which(df[, i] %in% max(df[, i]))
    commonFoundersRows = intersect(foundersLastRow, foundersThisRow)
    if (length(commonFoundersRows) > 0) {
      blocks = append(blocks, block)
    } else {
      # Confirm recombination events using the context
      # haplotype sequences
      
      # The context will be split into four regions: (1)
      # upwards, -1-k away from current haplotype, (2)
      # last haplotype, i-1, (3) current haplotype, i, (4)
      # downwards, +k away from current haplotype.
      # Comparing the upwards and the downwards context
      # regions will provide support or contradict the
      # putative crossover between the contrasting current
      # and next haplotypes
      upperLimit = ifelse(i-1-k >= 1, i-1-k, 1)
      bottomLimit = ifelse(i+k <= ncol(df), i+k, ncol(df))
      nextSNP = ifelse(i+1 <= ncol(df), i+1, ncol(df))
      upwardKwindow = as.matrix(df[, upperLimit:(i-2)])
      downwardKwindow = as.matrix(df[, nextSNP:bottomLimit])
      foundersUpwardKwindowRows = which(rowSums(upwardKwindow) %in% max(rowSums(upwardKwindow)))
      foundersDownwardKwindowRows = which(rowSums(downwardKwindow) %in% max(rowSums(downwardKwindow)))
      noHaplotypeSwitch = any(foundersUpwardKwindowRows %in% foundersDownwardKwindowRows)
      noHaplotypeSwitch = ifelse(k == 0, FALSE, noHaplotypeSwitch)
      
      # If any common founder is present both in the upwards
      # and downwards context regions, do not count as
      # block
      if (noHaplotypeSwitch) {
        blocks = append(blocks, block)
      } else {
        block = block + 1
        blocks = append(blocks, block)
        # Extract the genotypes within the start and end of
        # the putative CO interval detected
        candidateRegion = as.data.frame(
          t(genotypes[c(foundersLastRow[1], foundersThisRow[1]),
                      slidingWindows$binStart_snps[i-1]:slidingWindows$binEnd_snps[i]]))
        lastColNames = colnames(candidateRegion)
        colnames(candidateRegion) = c("upwardFounder", "downwardFounder")
        # Keep SNPs that are different between the two
        # recombining founders and then select the interval
        # between the first and the second
        candidateRegion$het = ifelse(candidateRegion$upwardFounder ==
                                       candidateRegion$downwardFounder,
                                     FALSE,
                                     TRUE)
        candidateRegion = candidateRegion %>% dplyr::filter(het)
        candidateRegion =
          candidateRegion %>%
          mutate(snp = rownames(candidateRegion),
                 hap = rleid(upwardFounder, downwardFounder))
        candidateRegion =
          suppressWarnings(suppressMessages(
            candidateRegion %>%
              group_by(hap) %>%
              slice(c(1, n())) %>%
              slice(c(2, 3)) %>%
              select(upwardFounder, snp, downwardFounder)
          ))
        fineInterval =
          paste0(lastColNames[1], ":", candidateRegion$snp[1], "-",
                 candidateRegion$snp[2], ":", lastColNames[2])
        fineIntervals = append(fineIntervals, fineInterval)
        fineBlocks = append(fineBlocks, block)
      }
    }
  }
  fineDF = data.frame(block = fineBlocks,
                      interval = fineIntervals)
  # Create a table with the haploblock coordinates and most
  # likely founders
  blocksDF = data.frame(block = blocks)
  blocksDF =
    suppressWarnings(suppressMessages(
      blocksDF %>%
        mutate(index = as.numeric(rownames(blocksDF))) %>%
        mutate(interval = colnames(df)[as.numeric(rownames(blocksDF))]) %>%
        group_by(block) %>%
        summarise(startIndex = first(index), endIndex = last(index),
                  startInterval = first(interval), endInterval = last(interval))
    ))
  blocksDF$startBP = as.numeric(str_match(blocksDF$startInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,3])
  blocksDF$endBP = as.numeric(str_match(blocksDF$endInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,2])
  blocksDF$startIntervalLength = as.numeric(str_match(blocksDF$startInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,3]) - as.numeric(str_match(blocksDF$startInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,2])
  blocksDF$endIntervalLength = as.numeric(str_match(blocksDF$endInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,3]) - as.numeric(str_match(blocksDF$endInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,2])
  # Generate a new haplotype path that assigns 1 to the most
  # likely founders within each block
  mostLikelyFoundersBlocks = c()
  mostLikelyFoundersBlocks2 = c()
  newDF = df[]
  blockSumScoreMax1 = c()
  blockSumScoreMax2 = c()
  blockSumScoreRest = c()
  for (i in 1:nrow(blocksDF)) {
    rowSumsBlock = rowSums(as.matrix(df[, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)]))
    bestFounderRows = which(rowSumsBlock %in% max(rowSumsBlock))
    bestFounderRows2 = which(rowSumsBlock %in% unique(sort(rowSumsBlock, decreasing = TRUE))[2])
    mostLikelyFounders = paste(rownames(df)[bestFounderRows], collapse = "/")
    mostLikelyFoundersBlocks = append(mostLikelyFoundersBlocks, mostLikelyFounders)
    mostLikelyFounders2 = paste(rownames(df)[bestFounderRows2], collapse = "/")
    mostLikelyFoundersBlocks2 = append(mostLikelyFoundersBlocks2, mostLikelyFounders2)
    newDF[, blocksDF$startIndex[i]:blocksDF$endIndex[i]] = 0
    newDF[bestFounderRows, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)] = 1
    blockSumScoreMax1 = append(blockSumScoreMax1, mean(rowSumsBlock[bestFounderRows])/(blocksDF$endIndex[i] - blocksDF$startIndex[i] + 1))
    blockSumScoreMax2 = append(blockSumScoreMax2, mean(rowSumsBlock[bestFounderRows2])/(blocksDF$endIndex[i] - blocksDF$startIndex[i] + 1))
    blockSumScoreRest = append(blockSumScoreRest, mean(rowSumsBlock[-c(bestFounderRows, bestFounderRows2)])/(blocksDF$endIndex[i] - blocksDF$startIndex[i] + 1))
  }
  blocksDF$mostLikelyFounders = mostLikelyFoundersBlocks
  blocksDF$mostLikelyFounders2 = mostLikelyFoundersBlocks2
  blocksDF$blockSumScore1 = blockSumScoreMax1
  blocksDF$blockSumScore2 = blockSumScoreMax2
  blocksDF$blockSumScoreRest = blockSumScoreRest
  mostLikelyFoundersExport = c()
  for (i in 1:ncol(newDF)) {
    mostLikelyFoundersExport = append(mostLikelyFoundersExport,
                                      paste(rownames(newDF)[which(as.matrix(newDF[, i]) %in% 1)], collapse = "/"))
  }
  exportList = list(as.data.frame(mostLikelyFoundersExport), fineDF, blocksDF)
  if (export == "mostLikelyFounders") {
    return (mostLikelyFoundersExport)
  } else if (export == "blocks") {
    return (blocksDF)
  } else if (export == "fineIntervals") {
    return (fineDF)
  } else if (export == "list") {
    return (exportList)
  }
}

#### inferMostLikelyPathAlternative5 ####
# The max-score founders on each block are compared btw the
# upwards and downwards windows and if there is any match, no
# crossover. This is an alternative more stringent approach
# to limit the huge number of crossovers obtained when the
# comparison was done with the highest-scoring founders of
# the whole windows
inferMostLikelyPathAlternative5 = function(df, k, slidingWindows, genotypes, export = "mostLikelyFounders") {
  # Tools for counting blocks
  blocks = c(1)
  block = 1
  # Tools for extracting putative CO intervals
  fineIntervals = c()
  fineBlocks = c()
  for (i in 2:ncol(df)) {
    # Compare the highest-score founders between the last
    # and the current position. If they are the same, keep
    # the same block, but if they differ, evaluate if the
    # recombination is true using the context
    foundersLastRow = which(df[, i-1] %in% max(df[, i-1]))
    foundersThisRow = which(df[, i] %in% max(df[, i]))
    commonFoundersRows = intersect(foundersLastRow, foundersThisRow)
    if (length(commonFoundersRows) > 0) {
      blocks = append(blocks, block)
    } else {
      # Confirm recombination events using the context
      # haplotype sequences
      
      # The context will be split into four regions: (1)
      # upwards, -1-k away from current haplotype, (2)
      # last haplotype, i-1, (3) current haplotype, i, (4)
      # downwards, +k away from current haplotype.
      # Comparing the upwards and the downwards context
      # regions will provide support or contradict the
      # putative crossover between the contrasting current
      # and next haplotypes
      upperLimit = ifelse(i-1-k >= 1, i-1-k, 1)
      bottomLimit = ifelse(i+k <= ncol(df), i+k, ncol(df))
      nextSNP = ifelse(i+1 <= ncol(df), i+1, ncol(df))
      upwardKwindow = as.matrix(df[, upperLimit:(i-2)])
      downwardKwindow = as.matrix(df[, nextSNP:bottomLimit])
      foundersUpwardKwindowRows = unlist(apply(upwardKwindow, 2, function(x) { which(x == max(x)) }))
      foundersDownwardKwindowRows = unlist(apply(downwardKwindow, 2, function(x) { which(x == max(x)) }))
      
      noHaplotypeSwitch = any(foundersUpwardKwindowRows %in% foundersDownwardKwindowRows)
      noHaplotypeSwitch = ifelse(k == 0, FALSE, noHaplotypeSwitch)
      
      # If any common founder is present both in the upwards
      # and downwards context regions, do not count as
      # block
      if (noHaplotypeSwitch) {
        blocks = append(blocks, block)
      } else {
        block = block + 1
        blocks = append(blocks, block)
        # Extract the genotypes within the start and end of
        # the putative CO interval detected
        candidateRegion = as.data.frame(
          t(genotypes[c(foundersLastRow[1], foundersThisRow[1]),
                      slidingWindows$binStart_snps[i-1]:slidingWindows$binEnd_snps[i]]))
        lastColNames = colnames(candidateRegion)
        colnames(candidateRegion) = c("upwardFounder", "downwardFounder")
        # Keep SNPs that are different between the two
        # recombining founders and then select the interval
        # between the first and the second
        candidateRegion$het = ifelse(candidateRegion$upwardFounder ==
                                       candidateRegion$downwardFounder,
                                     FALSE,
                                     TRUE)
        candidateRegion = candidateRegion %>% dplyr::filter(het)
        candidateRegion =
          candidateRegion %>%
          mutate(snp = rownames(candidateRegion),
                 hap = rleid(upwardFounder, downwardFounder))
        candidateRegion =
          suppressWarnings(suppressMessages(
            candidateRegion %>%
              group_by(hap) %>%
              slice(c(1, n())) %>%
              slice(c(2, 3)) %>%
              select(upwardFounder, snp, downwardFounder)
          ))
        fineInterval =
          paste0(lastColNames[1], ":", candidateRegion$snp[1], "-",
                 candidateRegion$snp[2], ":", lastColNames[2])
        fineIntervals = append(fineIntervals, fineInterval)
        fineBlocks = append(fineBlocks, block)
      }
    }
  }
  fineDF = data.frame(block = fineBlocks,
                      interval = fineIntervals)
  # Create a table with the haploblock coordinates and most
  # likely founders
  blocksDF = data.frame(block = blocks)
  blocksDF =
    suppressWarnings(suppressMessages(
      blocksDF %>%
        mutate(index = as.numeric(rownames(blocksDF))) %>%
        mutate(interval = colnames(df)[as.numeric(rownames(blocksDF))]) %>%
        group_by(block) %>%
        summarise(startIndex = first(index), endIndex = last(index),
                  startInterval = first(interval), endInterval = last(interval))
    ))
  blocksDF$startBP = as.numeric(str_match(blocksDF$startInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,3])
  blocksDF$endBP = as.numeric(str_match(blocksDF$endInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,2])
  blocksDF$startIntervalLength = as.numeric(str_match(blocksDF$startInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,3]) - as.numeric(str_match(blocksDF$startInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,2])
  blocksDF$endIntervalLength = as.numeric(str_match(blocksDF$endInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,3]) - as.numeric(str_match(blocksDF$endInterval, "chr\\w+_([0-9]+)-chr\\w+_([0-9]+)")[,2])
  # Generate a new haplotype path that assigns 1 to the most
  # likely founders within each block
  mostLikelyFoundersBlocks = c()
  mostLikelyFoundersBlocks2 = c()
  newDF = df[]
  blockSumScoreMax1 = c()
  blockSumScoreMax2 = c()
  blockSumScoreRest = c()
  for (i in 1:nrow(blocksDF)) {
    rowSumsBlock = rowSums(as.matrix(df[, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)]))
    bestFounderRows = which(rowSumsBlock %in% max(rowSumsBlock))
    bestFounderRows2 = which(rowSumsBlock %in% unique(sort(rowSumsBlock, decreasing = TRUE))[2])
    mostLikelyFounders = paste(rownames(df)[bestFounderRows], collapse = "/")
    mostLikelyFoundersBlocks = append(mostLikelyFoundersBlocks, mostLikelyFounders)
    mostLikelyFounders2 = paste(rownames(df)[bestFounderRows2], collapse = "/")
    mostLikelyFoundersBlocks2 = append(mostLikelyFoundersBlocks2, mostLikelyFounders2)
    newDF[, blocksDF$startIndex[i]:blocksDF$endIndex[i]] = 0
    newDF[bestFounderRows, seq(blocksDF$startIndex[i], blocksDF$endIndex[i], 1)] = 1
    blockSumScoreMax1 = append(blockSumScoreMax1, mean(rowSumsBlock[bestFounderRows])/(blocksDF$endIndex[i] - blocksDF$startIndex[i] + 1))
    blockSumScoreMax2 = append(blockSumScoreMax2, mean(rowSumsBlock[bestFounderRows2])/(blocksDF$endIndex[i] - blocksDF$startIndex[i] + 1))
    blockSumScoreRest = append(blockSumScoreRest, mean(rowSumsBlock[-c(bestFounderRows, bestFounderRows2)])/(blocksDF$endIndex[i] - blocksDF$startIndex[i] + 1))
  }
  blocksDF$mostLikelyFounders = mostLikelyFoundersBlocks
  blocksDF$mostLikelyFounders2 = mostLikelyFoundersBlocks2
  blocksDF$blockSumScore1 = blockSumScoreMax1
  blocksDF$blockSumScore2 = blockSumScoreMax2
  blocksDF$blockSumScoreRest = blockSumScoreRest
  mostLikelyFoundersExport = c()
  for (i in 1:ncol(newDF)) {
    mostLikelyFoundersExport = append(mostLikelyFoundersExport,
                                      paste(rownames(newDF)[which(as.matrix(newDF[, i]) %in% 1)], collapse = "/"))
  }
  exportList = list(as.data.frame(mostLikelyFoundersExport), fineDF, blocksDF)
  if (export == "mostLikelyFounders") {
    return (mostLikelyFoundersExport)
  } else if (export == "blocks") {
    return (blocksDF)
  } else if (export == "fineIntervals") {
    return (fineDF)
  } else if (export == "list") {
    return (exportList)
  }
}

#### RECORD CROSSOVER COORDINATES ####
# Record the first element of each haploblock as the crossover interval
getCrossovers = function(df) {
  df =
    df %>%
    mutate(interval = rownames(df), block = rleid(mostLikelyFounder)) %>%
    select(interval, mostLikelyFounder, block) %>%
    group_by(block) %>%
    slice(c(1,n()))
  df = df %>% dplyr::filter(mostLikelyFounder != "")
  return(df)
}

#### IMPUTE REPETITIVE HAPLOBLOCKS ####
# ABABACCC --> AAAAACCC
imputeRepetitiveSwitches = function(crossoversDF) {
  # crossoversDF = crossoversDF %>%
  #  dplyr::filter(!grepl("/", mostLikelyFounder),
  #                !grepl("/", mostLikelyFounder))
  founders = crossoversDF$mostLikelyFounder
  newFounders = founders
  for (i in seq(6, length(founders), 2)) {
    twoFoundersBack = newFounders[i-4]
    oneFounderBack = newFounders[i-2]
    thisFounder = newFounders[i]
    if (thisFounder == twoFoundersBack) {
      newFounders[i-2] = thisFounder
      newFounders[i-3] = thisFounder
    }
  }
  crossoversDF$mostLikelyFounder = newFounders
  crossoversDF =
    as.data.table(crossoversDF) %>%
    select(!block) %>%
    mutate(block = rleid(mostLikelyFounder))
  # Next loop substitutes block to merge haploblocks that
  # share founders (BP101 and BP101/BP103, for ex)
  # blocks = c(1)
  # block = 1
  # for (i in 2:nrow(crossoversDF)) {
  #   commonHaplotypes =
  #     str_extract(crossoversDF$mostLikelyFounder[i],
  #                 paste(strsplit(paste(crossoversDF$mostLikelyFounder[i-1], collapse = ", "), ", ")[[1]], collapse = "|"))
  #   condition = length(commonHaplotypes) > 0
  #   if (condition) {
  #     blocks = append(blocks, block)
  #     crossoversDF$mostLikelyFounder[i-1] = paste(commonHaplotypes, collapse = ", ")
  #     crossoversDF$mostLikelyFounder[i] = paste(commonHaplotypes, collapse = ", ")
  #   } else {
  #     block = block + 1
  #     blocks = append(blocks, block)
  #   }
  # }
  #crossoversDF$block = blocks
  crossoversDF =
    crossoversDF %>%
    group_by(block) %>%
    slice(c(1,n()))
  return(crossoversDF)
}

#### Print session info #####
#sessionInfo()
