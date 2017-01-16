#' Add Reactome hierarchical structure to rows of a matrix or data.frame of data for Reactome pathways
#'
#' @param sl_table Input data: matrix or data.frame with a row for each Reactome pathway. The first column ("Probe") must be the names of Reactome pathways and the second column ("Reactome ID") is used if given.
#' @param store_temp logical. Whether or not to store the pathway structure data. Defaults to FALSE (sourcing files from URL) unless files found in working directory.
#' @return The sum input data sorted by Reactome pathways with an additional column (left) for annotation.

add_tree_reactome <- function(sl_table, store_temp = (file.exists("ReactomePathwaysRelation.txt") & file.exists("ReactomePathways.txt"))){
  sl_table <- data.frame(sl_table, stringsAsFactors = F, row.names = 1:nrow(sl_table))

  ## Add organism info to pathway names
  reactomeName <- paste0("Homo sapiens: ",sl_table$Probe)

  ## Pathway names
  if(store_temp){
    ## download and import file
    if(!file.exists("ReactomePathways.txt")) curl_download("http://www.reactome.org/download/current/ReactomePathways.txt", "ReactomePathways.txt")
    pathway_info <- fread('ReactomePathways.txt', header=F, data.table = F)
  } else {
    ## import from web
    pathway_info <- fread('http://www.reactome.org/download/current/ReactomePathways.txt', header=FALSE, data.table = F)
  }
  ## Extract human data
  pathway_info <- pathway_info[grep("R-HSA-", pathway_info[,1]),]

  if("Reactome ID" %in%  colnames(sl_table) == T){
    #Retrieve known IDs
    reactomeID <- paste0("R-HSA-", sl_table$`Reactome ID`)
  } else {
    ## infer IDs from names
    ## Mapping from pathway names to short IDs
    pathway_ids <- as.list(reactomePATHNAME2ID)
    ## Convert pathway names to IDs, and add "R-HSA-" to each
    reactomeID <- paste0("R-HSA-",as.vector(unlist(lapply(pathway_ids[match(reactomeName, names(pathway_ids))], function(x) x[1]))))
  }

  if(exists("pathway_info")){
  ## What doesn't match between Tom's results, and the files from Reactome?
  ## (Haven't dealt with this issue)
  print(paste(sum(is.na(match(reactomeID, pathway_info[,1]))), "pathway names in dataset not found in reactome tree"))
  ##66 missing if IDs used
  print(paste(sum(is.na(match(pathway_info[,1], reactomeID,))), "reactome pathways not included in input dataset"))
  ##482 missing if names used
  }

  ## Reactome data:
  ## http://www.reactome.org/pages/download-data/
  ## Read in pathway relationship data file from reactome
  ## http://www.reactome.org/download/current/ReactomePathwaysRelation.txt
  if(store_temp){
    ## download and import file
    if(!file.exists("ReactomePathwaysRelation.txt")) curl_download("http://www.reactome.org/download/current/ReactomePathwaysRelation.txt", "ReactomePathwaysRelation.txt")
    pathway_str <- fread('ReactomePathwaysRelation.txt', header=FALSE, data.table = F)
  } else {
    ## import from web
    pathway_str <- fread('http://www.reactome.org/download/current/ReactomePathwaysRelation.txt', header=FALSE, data.table = F)
  }

  ## Find pathways in the up- and down-stream pathway relationship lists
  u <- match(reactomeID, pathway_str[,1])
  d <- match(reactomeID, pathway_str[,2])

  ## pathways not listed in hierarchy
  print(paste(table( reactomeID%in%pathway_str[,1] | (reactomeID%in%pathway_str[,2]) ), "pathways", c("not", ""), "listed"))

  ## Keep pathways that are reported in results
  kp <- pathway_str[,1]%in%reactomeID | pathway_str[,2]%in%reactomeID

  ## create format needed for generating tree
  pathway_str_kp <- pathway_str[kp,]
  names(pathway_str_kp) <- c("from","to")

  ## Get additional upstream pathways
  all_kp <- unique(c(as.vector(pathway_str_kp[,1]),as.vector(pathway_str_kp[,2])))
  top_kp <- all_kp[all_kp%in%pathway_str_kp[,1] & !(all_kp%in%pathway_str_kp[,2])]
  ## Create tree format data
  pathway_tree <- rbind(data.frame(from=rep("Reactome", length(top_kp)), to=top_kp), pathway_str_kp)

  ## Use data.tree package to generate tree
  treeDat <- FromDataFrameNetwork(pathway_tree)

  ## Convert to data.frame
  treeDf <- print(treeDat, limit=10000)

  ## Match data onto the tree
  matched_table <- list()
  for(i in 1:length(reactomeID)) matched_table[[i]] <- grep(paste0(reactomeID[i]," "),treeDf[,1])
  matched_names <- rep(" ", nrow(treeDf))
  for(i in 1:length(reactomeID)) matched_names[matched_table[[i]]] <- reactomeName[i]

  ## Get rid of some of the weird characters in the tree
  tree2 <- gsub("¦","",gsub("°","",gsub("--","",gsub(" ","",treeDf[,1]))))
  tree3 <- pathway_info[match(tree2,pathway_info[,1]),2]

  ## Build new data.frame with tree ordering
  ordered_tree <- data.frame(tree=treeDf[,1], path=tree3, tom=gsub("Homo sapiens: ","",matched_names))[-1,]
  ordered_treeDF <- as.data.frame(matrix(NA, nrow(ordered_tree), ncol(sl_table)))
  for(i in 1:length(reactomeID)){
    if(length(matched_table[[i]]>0)){
      for(j in 1:length(matched_table[[i]])) ordered_treeDF[matched_table[[i]][j],] <- sl_table[i,]
    }
  }
  colnames(ordered_treeDF) <- colnames(sl_table)
  full_tree <- cbind(ordered_tree,ordered_treeDF[-1,])
  full_tree[,1] <- gsub(" ¦"," ¦",full_tree[,1])
  crop_tree <- full_tree[,-c(3:4)]
  return(crop_tree)
}
