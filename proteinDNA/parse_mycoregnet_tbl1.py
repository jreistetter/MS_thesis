# Script to remove the gene names from TF targets listed in
# table 1 of MycoRegNet article. Save them as a tab-delimited file
# where first column is regulator, second column is comma-delimited
# list of its targets.
# 
# Written by Joe Reistetter

#laptop path
root_path = "/Users/jreistetter/schoolDB/Dropbox/thesis_work/data/protein-DNA/MycoRegNet/"
#OHSU path
#root_path = "~/Dropbox/thesis_work"

path = root_path+"mycoregnet_tbl1_hand_cleaned.txt"
lines = open(path, 'r').readlines()
out_f = open(root_path+"mycoregnet_tbl1_parsed.txt", 'w')

clean = []
for line in lines:
    line = line.split('\t')
    reg = line[0]
    targets = ','.join([chunk.strip().split(' ')[0] for chunk in line[1].split(',')])
    out_f.write('\t'.join([reg,targets]) + '\n')

out_f.close()