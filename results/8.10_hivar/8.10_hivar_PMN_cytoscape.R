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

# Load pathways to create nodes for pathways from parent to target TF
setwd("~/Dropbox/thesis_work/PMN_output/")
paths.raw <- read.table(
  "8.10.12_30_mods_0.1_highvar/8.10.12_30_mods_0.1_highvar_mods_pathways.txt",
  head=F, sep=":")
colnames(paths.raw) <- c("moduleID", "pathway")
dim(paths.raw)

# Select the pathways for each mod. Write an edge from parent to the module
# it belongs to. Then write an edge from parent to first protein in pathway,

setwd("~/Dropbox/thesis_work/data/8.10_highvar_results/cytoscape/")
out_f <- file("PMN_pathways.sif", "w")
write("source\tinteraction\ttarget", out_f)

for (ID in modIDs){
  # Get all pathways associated with module
  mod.paths <- paths.raw[paths.raw$moduleID==ID,]
  
  for (i in 1:nrow(mod.paths)){
    # Loop through each pathway
    path <- unlist(strsplit(mod.paths[i,2], "\t", fixed=T))
    # Write out edge from module (that parent belongs to) to parent
    parent.mod <- gene.mod(path[1], modules)
    if (length(parent.mod) > 0){
      line <- paste(parent.mod, "mp", path[1], sep="\t")
      write(line, out_f)
    }
    path.length <- length(path)
    # Write out edges for pathway
    for (i in 1:(path.length-1)){
      line <- paste(path[i], "pp", path[i+1], sep="\t")
      write(line, out_f)
    }
    
    # Write out last edge from TF to parents child module
    line <- paste(path[path.length], "pm", ID, sep="\t")
    write(line, out_f)
  }
}

close(out_f)

## Write out the pathways again, but with no links to the modules

setwd("~/Dropbox/thesis_work/data/8.10_highvar_results/cytoscape/")
out_f <- file("PMN_pathways_PPI_only.sif", "w")
write("source\tinteraction\ttarget", out_f)

for (ID in modIDs){
  # Get all pathways associated with module
  mod.paths <- paths.raw[paths.raw$moduleID==ID,]
  
  for (i in 1:nrow(mod.paths)){
    # Loop through each pathway
    path <- unlist(strsplit(mod.paths[i,2], "\t", fixed=T))
    path.length <- length(path)
    # Write out edges for pathway
    for (i in 1:(path.length-1)){
      line <- paste(path[i], "pp", path[i+1], sep="\t")
      write(line, out_f)
    }

  }
}

close(out_f)

# Write out a node attribute file with module or protein attribute

# Read the network edges back in to get list of all nodes
net <- read.table("PMN_pathways.sif", head=T, sep="\t")

all_nodes <- unique(c(net$source, net$target))
length(all_nodes)
# [1] 102

mod.nodes <- all_nodes[grepl("mod", all_nodes, fixed=T)]
length(mod.nodes)
# [1] 28

gene.nodes <- all_nodes[grepl("RV", all_nodes, fixed=T)]
length(gene.nodes)
# [1] 73

node.attributes <- data.frame(nodeID=c(mod.nodes, gene.nodes))
node.attributes$type <- "module"
node.attributes[grepl("RV", node.attributes$nodeID),2] <- "protein"

nrow(node.attributes)
# [1] 101 = 28 + 73

write.table(node.attributes, "PMN_node_attributes.noa",
            col.names=T, row.names=F, quote=F, sep="\t")











