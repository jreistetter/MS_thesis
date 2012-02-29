# Script to create the protein-protein edge list for input to learn a PMN. 
# 
# 1. Import protein-protein interactions from two datasets:
#   -STRING Mtb H37Rv functional protein associations from flatfile
#   -Protein-protein interactions derived from bacterial 2-hybrid assays by Wang et al. 2010 PMID: 20973567
# 
# 2. Filter STRING associations for confidence > 0.9
# 
# 3. Calculate the confidence score for each edge as given in Novershtern et al. 2011 PMID: 21685068 
# 
# 4. Write the edges and their confidence scores to a file readable by the PMN software from http://www.compbio.cs.huji.ac.il/PMN/
# format:
# tab-delimited with 4 columns <prot1> <prot2> <p-value> <direction>
# direction = 0 (undirected), 1 (directed forward), or 2 (directed backwards) 
# In pp interaction files, it will typically be 0 (undirected). 

