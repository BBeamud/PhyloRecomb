library(sqldf)
library("optparse")

option_list = list(
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="Sample name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
pattern1=paste("*", opt$sample, sep='')
pattern2=paste(pattern1, "*.summary", sep='')
command = paste("ls ", pattern2, sep='')
# Get all results for this sample
summs <- try(system(command, intern = TRUE))
#summs = list.files(pattern = pattern2)
# Create empty dataframe
matriz=matrix(data= NA, length(summs), 3)
final = as.data.frame(matriz)


for (i in 1:length(summs)){
  table = read.table(summs[i], header=T)
  # We select ELW and AU test plus subtype and fragment information
  main = table[,c(11, 13, 14, 15)]
  # We create BED file proper to draw CRF in Los Alamos Drawing Tool
  interval =strsplit(as.character(main$Fragment), "-")[[1]]
  final[i,"V1"] = as.numeric(interval[[1]])
  final[i,"V2"] = as.numeric(interval[[2]])
  pos = rowSums(main == "+")
  # Both tests has to be positive to consider verosimil
  ## AJUSTAR LUEGO MEJORRR CONSIDERANDO DIFERENCIAS ENTRE VEROSIMILITUDES!!
  index = which(pos == 2)
  # If only one subtype is verosimil
  if (length(index) == 1) {
    final[i,"V3"] = as.character(main$Subtype[[index]])
  }
  # If more than one subtype is verosimil we have to include all
  else {
    subtypes= sapply(subtypes, function(x){
      # we had dto add '/' between subtypes except in last one
      if (x != subtypes[length(subtypes)]) {
        subtypes[x]=paste(x,'/',sep='')
      }
      else {
        subtypes[x]=x
      }
    })
    # we collapse the list in one string to save it 
    subtypes=paste(subtypes, collapse = '')
    final[i,"V3"]=subtypes
    }
  }


# Sort data frame by start position
final=final[order(final$V1),]
output=paste(opt$sample, ".breakpoints", sep='')
write.table(final, output, col.names = F, row.names = F, quote = F)





