# source("https://bioconductor.org/biocLite.R") 
# biocLite('org.Hs.eg.db')
geneIDannotation <- function(geneLists=c(1,2,9),name=T,map=T,ensemble=F,accnum=F){
  ## input ID type : So far I just accept entrezID or symbol
  ## default, we will annotate the entrezID and symbol, chromosone location and gene name 
  
  suppressMessages(library("org.Hs.eg.db"))
  all_EG=mappedkeys(org.Hs.egSYMBOL) 
  EG2Symbol=toTable(org.Hs.egSYMBOL)
  if( all(! geneLists %in% all_EG) ){
    inputType='symbol'
    geneLists=data.frame(symbol=geneLists)
    results=merge(geneLists,EG2Symbol,by='symbol',all.x=T)
  }else{
    inputType='entrezID'
    geneLists=data.frame(gene_id=geneLists)
    results=merge(geneLists,EG2Symbol,by='gene_id',all.x=T)
  }
    
  if ( name ){
    EG2name=toTable(org.Hs.egGENENAME)
    results=merge(results,EG2name,by='gene_id',all.x=T)
  }
  if(map){
    EG2MAP=toTable(org.Hs.egMAP)
    results=merge(results,EG2MAP,by='gene_id',all.x=T)
  }
  if(ensemble){
    EG2ENSEMBL=toTable(org.Hs.egENSEMBL)
    results=merge(results,EG2ENSEMBL,by='gene_id',all.x=T)
  }
  if(accnum){
    EG2accnum=toTable(org.Hs.egREFSEQ) 
    results=merge(results,EG2MAP,by='gene_id',all.x=T)
  }
  return(results)
}
b = geneIDannotation(Module_Gene[,2])
b = b[,1:2]
b = unique(b)

# ID_Change=geneIDannotation(Symbol_list)