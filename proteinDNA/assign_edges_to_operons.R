# Script to take regulator-target lists and operon data from tbdb.org
# and parse it into a list of protein-DNA interactions for use in the PMN software.
# 
# The regulator-target data has the start and stop coords of the peak
# 
# The operon data has:
#   the start and stop coords of the operon
#   the strand of the operon
#   the operon name, which isn't very useful for determining membership
#   the length
# 
# Written by Joe Reistetter

# -import both data
# -get H37Rv annotation from biomaRt
# -use gene start/stop and operon start/stop coords to assign genes to operons
# -for each gene in an operon, create an interaction between regulator and gene (protein-DNA edge)