#!/usr/bin/env Rscript
# This program create a new alternative tree for each alternative subtype.
# No bootstrap accepted. Assumes that reference subtype is the topology of that fragment.  
# Better without distances because alternative trees doesnt have
# Usage: Rscript --vanilla alternative_trees.R -t name.tre -s tip_name -a tip_name,tip_name,... -p prefix

library(ape)
library(phytools)
library("optparse")

option_list = list(
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="Reference tree name", metavar="character"),
  make_option(c("-r", "--ref"), type="character", default=NULL, 
              help="Assigned subtype", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="Name of recombinant sample", metavar="character"),
  make_option(c("-a", "--alt"), type="character", default=NULL, 
              help="List of alternative subtypes", metavar="character"),
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="Prefix of the output file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# get reference tree 
#tree=read.tree("79_fff_1920-3154_subaln.tre")
# get alternatives subtypes
#subtypes = "C,F1"
#ref = 'B'
# get name prefix for output
#name="pepe"
# First remove sample
#r_sample = '79_fff_1920-3154'
#int_tree = drop.tip(tree, r_sample)
# Save trees



# get reference tree 
tree=read.tree(opt$tree)
ref=opt$ref
# get alternatives subtypes
subtypes=opt$alt
subtypes=as.list(strsplit(subtypes, ",")[[1]])
# get name prefix for output
name=opt$prefix
output = paste(name, ".trees", sep='')


# Sample to remove
r_sample=opt$sample
int_tree = drop.tip(tree, r_sample)
# write ref tree according with assigned subtype
node_ref = which(int_tree$tip.label==ref)
ref_tree = bind.tip(int_tree, tip.label=r_sample, where=node_ref)
write.tree(ref_tree, output)



lapply(subtypes, function(x){
  # Get nodes from alternative subtypes
  node_alt = which(int_tree$tip.label==x)
  # Add sample to that node
  alt_tree = bind.tip(int_tree, tip.label=r_sample, where=node_alt)
  # Add alt tree to all tree file 
  write.tree(alt_tree, output, append = TRUE)
})


