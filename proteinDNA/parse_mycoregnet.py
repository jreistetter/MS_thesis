# Script to parse regulator-target pairs from
# tables in MycoRegNet article. Save them as a tab-delimited file
# where first column is regulator, second column is target
# 
# Written by Joe Reistetter

#laptop path
root_path = "/Users/jreistetter/schoolDB/Dropbox/thesis_work/data/protein-DNA/MycoRegNet/"
#OHSU path
root_path = "/Domain/ohsum01.ohsu.edu/Users/reistett/Dropbox/thesis_work/data/protein-DNA/MycoRegNet/"

#Table 1
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

#Table 2
path = root_path+"mycoregnet_tbl2_hand_cleaned.txt"
#Read in lines, skip first header row
lines = open(path, 'rU').readlines()[1:]
out_f = open(root_path+"mycoregnet_tbl2_parsed.txt", 'w')

for line in lines:
	line = line.strip().split('\t')
	regulator = line[0]
	
	#Some targets are part of an operon and have a third column listing
	#the operon members, delimited by '-'.
	if len(line) == 3:
		targets = line[2].split('-')
		
		for target in targets:
			out_f.write('\t'.join([regulator, target])+'\n')
	
	#If no operon associated, then it is just a regulator/target pair
	else:
		target = line[1]
		out_f.write('\t'.join([regulator, target])+'\n')

out_f.close()

			