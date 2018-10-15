# Resume results 
library("optparse")
option_list = list(
  make_option(c("-s", "--summary"), type="character", default=NULL, 
              help="Summary of results", metavar="character"),
  make_option(c("-r", "--ref"), type="character", default=NULL, 
              help="Assigned subtype", metavar="character"),
  make_option(c("-a", "--alt"), type="character", default=NULL, 
              help="List of alternative subtypes", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
results = read.table(opt$summary, col.names = c("ELW", "AU"), header=T)

ref_sub=opt$ref
alt_sub=opt$alt
subtypes=as.list(strsplit(alt_sub, ",")[[1]])
all_s=as.vector(c(ref_sub,subtypes))
results["Subtype"]=unlist(all_s)
# Based only in ELW
accept = results[which(as.numeric(results$ELW) >0.05),]
accept = accept[order(accept$ELW, decreasing = T),]
best=accept$Subtype[1]
if (nrow(accept)<2){
  print(accept$Subtype)
} else {
  all = best
  for (i in 2:nrow(accept)){
    dif = accept$ELW[1] - accept$ELW[i]
    if (dif > 0.5){
      subtype = tolower(accept$Subtype[i])
      all = c(all, subtype)
        } else {
          subtype = accept$Subtype[i]
          all = c(all, subtype)
        }
      }
  } 
  print(paste(all, collapse="|"))
  
