### Lada's function

#source("https://bioconductor.org/biocLite.R")
#biocLite("GO.db")

library(GO.db)

obj.exist <- function(object){
  exists(as.character(substitute(object)))
}


reduce_GO = function(godf,ngene_lim = 500){ 
  ### godf = dataframe with at least these column names "Geneset.ID" and "N.genes" - no particular order necessary, other columns ok so long as at least have the required ones
  ### ngene_lim is to pre-filter gene sets so that only smaller, more informative gene sets are used
  ### VERY IMPORTANT that the input godf is already ordered by some value, like FDR, rank, etc
  ### I recommend doing reduce_GO twice if results from your gene set enrichment test is 2-sided (e.g. enriched/depleted, up/down)
  
  godf = subset(godf, N.genes <=ngene_lim)
  
  if (!obj.exist(goall.off) & !obj.exist(goall.par) & !obj.exist(goall.anc)){
    # find related GO terms, look at offspring and at parents
    library(GO.db)
    gobp.off <- as.list(GOBPOFFSPRING)
    # Remove GO IDs that do not have any offspring
    gobp.off <- gobp.off[!is.na(gobp.off)]
    
    gobp.par <- as.list(GOBPPARENTS)
    # Remove GO IDs that do not have any parents
    gobp.par <- gobp.par[!is.na(gobp.par)]
    
    gobp.anc <- as.list(GOBPANCESTOR)
    # Remove GO IDs that do not have any ancestors
    gobp.anc <- gobp.anc[!is.na(gobp.anc)]
    
    
    gocc.off <- as.list(GOCCOFFSPRING)
    # Remove GO IDs that do not have any offspring
    gocc.off <- gocc.off[!is.na(gocc.off)]
    
    gocc.par <- as.list(GOCCPARENTS)
    # Remove GO IDs that do not have any parents
    gocc.par <- gocc.par[!is.na(gocc.par)]
    
    gocc.anc <- as.list(GOCCANCESTOR)
    # Remove GO IDs that do not have any ancestors
    gocc.anc <- gocc.anc[!is.na(gocc.anc)]
    
    
    gomf.off <- as.list(GOMFOFFSPRING)
    # Remove GO IDs that do not have any offspring
    gomf.off <- gomf.off[!is.na(gomf.off)]
    
    gomf.par <- as.list(GOMFPARENTS)
    # Remove GO IDs that do not have any parents
    gomf.par <- gomf.par[!is.na(gomf.par)]
    
    gomf.anc <- as.list(GOMFANCESTOR)
    # Remove GO IDs that do not have any ancestors
    gomf.anc <- gomf.anc[!is.na(gomf.anc)]
    
    goall.par = c(gobp.par,gocc.par,gomf.par)
    goall.off = c(gobp.off,gocc.off,gomf.off)
    goall.anc = c(gobp.anc,gocc.anc,gomf.anc)
    
  }
  
  if( ! "Description" %in% names(godf)){
    godf$Description = Term(godf$Geneset.ID)
  }
  golist = godf$Geneset.ID
  godf = data.frame(orig.row.num = 1:nrow(godf),godf)
  ranked.list = godf[1,]
  for (i in 2:length(golist)){
    gterm = golist[i]
    related = c(unlist(goall.off[which(names(goall.off)==gterm)]),unlist(goall.par[which(names(goall.par)==gterm)]),unlist(goall.anc[which(names(goall.anc)==gterm)]))
    
    uniq = ifelse(sum(golist[1:(i-1)] %in% related)>0,0,1)
    if(uniq==1){
      ranked.list = rbind(ranked.list,godf[i,])
    }
  }
  
  return(unique(ranked.list[complete.cases(ranked.list),]))
}
