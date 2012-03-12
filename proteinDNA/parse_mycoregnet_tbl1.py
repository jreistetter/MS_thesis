# Script to parse regulator-target pairs from
# table 1 of MycoRegNet article. Save them as a tab-delimited file
# where first column is regulator, second column is target
# 
# Written by Joe Reistetter

#laptop path
root_path = "/Users/jreistetter/schoolDB/Dropbox/thesis_work/data/protein-DNA/MycoRegNet/"
#OHSU path
#root_path = "~/Dropbox/thesis_work"

path = root_path+"mycoregnet_tbl1_hand_cleaned.txt"
lines = open(path, 'r').readlines()
out_f = open(root_path+"mycoregnet_tbl1_parsed.txt", 'w')

for line in lines:
    line = line.split('\t')
    reg = line[0]
    targets = [chunk.strip().split(' ')[0] for chunk in line[1].split(',')]
    for target in targets:
        out_f.write('\t'.join([reg,target]) + '\n')

out_f.close()