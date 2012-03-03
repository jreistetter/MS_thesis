# Script to take expression data stored as RData objects, discretize them, and then write the data to a file readable by the PMN software.
# 
# 1 - Load expression data from RData objects
# 2 - Discretize data at a given fold change
# 3 - Merge the data frames from each experiment
# 4 - Write the combined dataframe to a PMN file