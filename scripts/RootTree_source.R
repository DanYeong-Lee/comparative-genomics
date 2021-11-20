# Title     : Root or unroot
# Objective : newick tree
# Created by: ldy93
# Created on: 2020-10-24


library(ape)
tree <- read.tree('./SpeciesTree/no_branch_length.nwk')

if(is.rooted(tree) == TRUE){
  unrooted_tree <- unroot(tree)
  write.tree(unrooted_tree, './SpeciesTree/unrooted.nwk')
  write.tree(tree, './SpeciesTree/rooted.nwk')
}else if(is.rooted(tree) == FALSE){
  rooted_tree <- root(tree, outgroup = '>>>>', resolve.root = TRUE)
  write.tree(rooted_tree, './SpeciesTree/rooted.nwk')
  write.tree(tree, './SpeciesTree/unrooted.nwk')
}


