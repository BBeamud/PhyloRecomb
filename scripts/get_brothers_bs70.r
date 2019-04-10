# Script to get closest taxa to one sample given a phylogeny 
library(phytools)
library(ape)
library(optparse)

option_list = list(
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="Reference tree name", metavar="character"),
  make_option(c("-s", "--sample"), type="character", default=NULL, 
              help="Name of recombinant sample", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
sample = opt$sample
tree_name = opt$tree 
tree=read.tree(tree_name)





# tree$edge is a matrix. First column is internal node and second external tips or other internal
# nodes connected with the formers.
# We get the tip number of our sample of interest 
tip_n=which(tree$tip.label==sample)
internal_node = getParent(tree, tip_n)
bs = as.numeric(tree$node.label[(internal_node-1)-(tree$Nnode+1)])
initial_n = tip_n
# We get the position in the matrix
#tip_index= which(tree$edge[,2]==tip_n)
# We get internal node
#internal_node = tree$edge[tip_index,1]
# We need to filter by bootstrap; to get the most recent clade with at least support >= 95 %. 
while (bs < 0.7) {
  internal_node = getParent(tree, tip_n)
  bs = as.numeric(tree$node.label[(internal_node-1)-(tree$Nnode+1)])
  tip_n = internal_node
}



# This function returns the set of nodes & tip numbers descended from closest node to our sample
descendants = getDescendants(tree, internal_node)
# We want only external tips, those that are not in the first column of the edge matrix  
descendants_tips = descendants[which(!descendants %in% tree$edge[,1])]
# Tips that are not my sample
descendants_tips= descendants_tips[descendants_tips != initial_n]
# We save tip names of all descendants 
taxas= character(0)
for (i in 1:length(descendants_tips)){
  name =  tree$tip.label[descendants_tips[i]]
  taxas = c(taxas,name)
}

# We can save all closest relative insolates, or the one with less distance to our sample
# A patristic distance is the sum of the lengths of the branches that link two nodes in a tree

for (i in 1:length(taxas)){
  if (i==1){
    dist0 = fastDist(tree, sample, taxas[i])
  } else {
    dist = fastDist(tree, sample, taxas[i])
    if (dist == dist0){
      dist0 = c(dist0, dist)
      closer = c(closer, taxas[i])
    }
    if (dist < dist0){
      dist0 = dist
      closer = taxas[i]
    # If both taxa present the same distance, we only retrieved first
  }
}
}
print(rev(taxas)) # We order the taxa; first closer. 

# Bootstrap is stored as node labels, they are actually associated with edges (i.e., splits) not with nodes
# Get bootstrap given an internal node:
#tree$node.label[(internal_node-1)-(tree$Nnode+1)]


