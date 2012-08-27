# Export the modules to cytoscape.
# The cytoscape network will show modules as nodes, parents as separate nodes,
# and indicate if a parent belongs to another module.
#
# TODO: add edges of PPI from parent to its child module

options(stringsAsFactors=F)
source("~/Dropbox/thesis_work/code/thesis_funcs.R")

# Load in module membership and module parents

setwd("~/Dropbox/thesis_work/PMN_output/8.10.12_30_mods_0.1_highvar/")

modules <- read.table("output/8.10_hivar_modules.txt")
modIDs <- unique(modules$moduleID)

cpds <- read.table("output/8.10.12_30_mods_0.1_highvar_CPDs.txt",
                   head=T, sep="\t")


# Check the parents of each module to see if they are in another module.
# If so, write out in cytoscape format: source | interaction | target | ..

setwd("~/Dropbox/thesis_work/data/8.10_highvar_results/cytoscape/")
out_f <- file("PMN_parent_module_membership.sif", "w")
write("source\tinteraction\ttarget\tparent_gene", out_f)

# Loop through all modules. Get the module parents, and for each parent
# check if it is in another module. If so, write out an edge from the
# module where the parent is a member to child module of the parent.
for (ID in modIDs){
  mod.parents <- get_parents(ID, cpds)
  
  for (parent in mod.parents){
    #Try to get the module ID for the parent
    parent.module <- modules[modules$gene==parent,]$moduleID
    
    #Write out edge if module found
    if (length(parent.module) > 0){
      line <- paste(parent.module, "mp", ID, parent, sep="\t")
      write(line, out_f)
    }
  }
}

close(out_f)