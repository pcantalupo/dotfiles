search_for_function = function (query, lib) {
  package = paste0("package:", lib)
  pacman::p_load(char = lib)
  grep(query, ls(package), value = TRUE, ignore.case = TRUE)
}




# Testing
# 3 human genes and convert to mouse
#myresults = getHomologousSymbols()
# 2 mouse genes and convert to human (Pbsn does not have a human homolog)
#myresults = getHomologousSymbols(symbols = c("Trp53", "Pbsn"), current = "mouse", target = "human")

# Downstream usage
# If you want to only keep the rows in the Results table that have a Target species symbol, do the following
#keep = !is.na(myresults$results[,2]) & myresults$results[,2] != ""
#myresults$results[keep,]

# Rat species is not yet supported
# Use 'force=TRUE' if you want to force connection to Ensembl Mart even if you already have a saved Mart.rds file for the current species
getHomologousSymbols = function(symbols = c("TP53", "RB1", "FOXP3"), current = "human", target = "mouse", force = FALSE) {
  # inspiration for this function comes from Aaron Lun: https://support.bioconductor.org/p/9136905/
  require('biomaRt')

  argg <- c(as.list(environment()))  # output function paramaters
  print(argg)
  
  toReturn = list()   # list to hold all the important objects created below
  
  # Ensembl datasets
  # First we need to pick an Ensembl dataset to use. You can see the list of datasets using the code below
      #ensemblmart <- useMart("ensembl")
      #datasets = listDatasets(ensemblmart)
      #head(datasets)
      #datasets[grep("mouse", datasets$description, ignore.case = TRUE),]
      #datasets[grep("rat", datasets$description, ignore.case = TRUE),]

  # The species datasets were determined using the above code
  species_datasets = c(human = "hsapiens_gene_ensembl", mouse = "mmusculus_gene_ensembl", rat = "rnorvegicus_gene_ensembl")
  current_dataset = as.character(species_datasets[names(species_datasets) == current])
  message(current_dataset)

  # Create (or load existing) ENSEMBL Mart
  ensemblmartrds = paste0("ensemblmart_", current, ".rds")
  if (file.exists(ensemblmartrds) && !force) {
    message("Loading ensemblmart from RDS file")
    ensemblmart = readRDS(ensemblmartrds)
  } else {
    message("Creating ensemblmart with 'useMart()'")
    # using uswest because of https://support.bioconductor.org/p/9144682/  # host default is www.ensembl.org
    ensemblmart = useMart(biomart = "ensembl", dataset = current_dataset, 
                          host="https://uswest.ensembl.org", verbose = TRUE)
  }
  ensemblmart
  saveRDS(ensemblmart, ensemblmartrds)

  # Attributes - features
  # Next, I used the listAttributes code below to determine the attribute values for the organism-specific symbols. At this step, we need to get the corresponding ENSEMBLID for the current species Symbols. We need the ENSEMBLIDs since they allow the connection to homolog information. At this step however, we cannot get homolog information yet because this information is on a different Attribute 'Page' and requires an independent getBM() function call.
      # ensembl_gene_id = EnsemblID (i.e. ENSGXXXXX)
      # hgnc_symbol = Human symbol
      # mgi_symbol = Mouse symbol
      # rgd_symbol = Rat symbol
      #grep("symbol", listAttributes(ensemblmart)[,1], value=TRUE) # show symbol attributes

  # The species symbol attributes were determined using the above code
  species_symbol_attributes = c(human = "hgnc_symbol", mouse = "mgi_symbol", rat = "rgd_symbol")
  current_symbol_attr = as.character(species_symbol_attributes[names(species_symbol_attributes) == current])
  attrs = c(current_symbol_attr, "ensembl_gene_id")
  message("\nRunning getBM to obtain ENSEMBLID for the Symbols (attributes are: ", attrs, ")")
  symbol2ensemblid <- getBM(attributes = attrs, filters = current_symbol_attr, values = symbols, mart = ensemblmart)
  toReturn[['symbol2ensemblid']] = symbol2ensemblid

  # Attributes - homologs
  # Third, using the current species ENSEMBLID, we can get the homologous ENSEMBLID and Symbol for the target species. I used the listAttributes code below to determine the attribute values for the organism-specific homologous information 
      # hsapiens_homolog_ensembl_gene = Human homologous EnsemblID (i..e ENSGXXXX)
      # hsapiens_homolog_associated_gene_name = Human homologous Symbol (i.e. TP53)
      # mmusculus_homolog_ensembl_gene = Mouse homologous EnsemblID (i..e ENSGMXXXX)
      # mmusculus_homolog_associated_gene_name = Mouse homologous Symbol (i.e. Trp53)
      # rnorvegicus_homolog_ensembl_gene = Rat homologous ENSEMBLID
      # rnorvegicus_homolog_associated_gene_name = Rat homologous Symbol

  # The species homolog attributes were determined using the above code
  species_homolog_attributes = list(human = c("hsapiens_homolog_ensembl_gene",
                                              "hsapiens_homolog_associated_gene_name"),
                                    mouse = c("mmusculus_homolog_ensembl_gene",
                                              "mmusculus_homolog_associated_gene_name"))
  current_homolog_attr = species_homolog_attributes[names(species_homolog_attributes) == target][[1]]
  attrs = c("ensembl_gene_id", current_homolog_attr)
  message("\nRunning getBM to obtain homologous ENSEMBLID and Symbols for target species (attributes are: ", attrs, ")")
  # in case of duplicate ENSEMBLIDs for the current species, we only use the unique values
  homolog_mapping = getBM(attributes = attrs,
                          filters="ensembl_gene_id", values=unique(symbol2ensemblid$ensembl_gene_id),
                          mart=ensemblmart)
  toReturn[['homolog_mapping']] = homolog_mapping
  
  # Merging getBM results
  # Fourth, we merge the two getBM() result tables. There can be duplicate current species ENSEMBLIDs since it may map to multiple target ENSEMBLIDs. Target ENSEMBLIDs may or may not have an associated gene symbol
  message("\nMerging the symbol2ensemblid and homolog_mapping tables")
  mapping = merge(symbol2ensemblid, homolog_mapping)
  toReturn[['mapping']] = mapping
  
  # Create results table with Symbols from Current species and Target species
  # Fifth, we create a dataframe in the same order as the current species Symbols and add the homologous target species symbols in the 2nd column. The tricky conceptual part is the 'match()' function that makes sure the 'mapping' table (from the 'merge' step above) is ordered the same as the current species Symbol vector 
  message("\nCreating final results table with Current and Target species Symbol values")
  current_homolog_symbol_attr = current_homolog_attr[2]
  targetsymbols = mapping[match(symbols, mapping[,current_symbol_attr]), current_homolog_symbol_attr]
  results = data.frame(symbols, targetsymbols)
  colnames(results) = c(current, target)
  toReturn[['results']] = results
  
  return(toReturn)
}


myRinfo = function() {

  message("\nIs Bioconductor valid?")
  print(suppressMessages(BiocManager::valid()))
  
  message("\nBioconductor version:")
  print(BiocManager::version())
  
  message("\nLibrary paths .libPaths():")
  print(.libPaths())
  
  message("\nSeurat version in .libpaths():")
  print(utils::packageVersion('Seurat'))
  
  locallib = "~/projects/Seurat/packages"
  message("\nChecking for Seurat in local library: ", locallib)
  if (dir.exists(locallib)) {
    res = utils::packageDescription("Seurat", lib.loc = locallib)
    print(paste0("Seurat version local: ", res[["Version"]]))
  } else {
    print(paste0("Local library ", locallib, "does not exist"))
  }
  
  message("\nVersion of select packages")
  pkgs = utils::installed.packages()
  p = c("SingleR", "celldex", "Seurat", "monocle3", "scater", "scran","Azimuth")
  print(pkgs[rownames(pkgs) %in% p,"Version", drop = FALSE])
  
  message("\nR version:")
  print(R.version.string)

}
