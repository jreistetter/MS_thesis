# Script to parse the screen scrape from operonDB.
# OperonDB doesn't have a function to download operons,
# but you can just copy and paste the screen into a text file.
# The data is tab delimited, with a few blank lines and descriptive
# lines that can be filtered out.
#
#
# Data was obtained from OperonDB:
# (http://operondb.cbcb.umd.edu/cgi-bin/operondb/operons.cgi)
# from the H37Rv page:
# http://operondb.cbcb.umd.edu/cgi-bin/operondb/pairs.cgi?genome_id=613
# Page was loaded, text was highlighted, copied, and pasted into a text file
# and saved as mtb_scrape.txt
#
# Written by Joe Reistetter

#laptop WD

operonDB_path = """/Users/jreistetter/schoolDB/Dropbox/thesis_work\
/data/protein-DNA/operonDB/mtb_scrape.txt"""

save_path = """/Users/jreistetter/schoolDB/Dropbox/thesis_work\
/data/protein-DNA/operonDB/operonDB_operons.txt"""

lines = open(operonDB_path, 'r').readlines()

goodlines = filter(lambda x: x.find('Rv') > -1, lines)

out_f = open(save_path, 'w')
header = '\t'.join(['rv1', 'rv2', 'conf', 'n_genomes'])
out_f.write(header + '\n')


for line in goodlines:
    #split into 3 tab-delimited columns
    sep = line.strip().split('\t')
    rv1 = sep[0].split()[0]
    rv2 = sep[1].split()[0]
    conf = sep[2].split()[0].split('=')[1]
    n_genomes = sep[2].split()[1].split('=')[1]
    
    line_out = '\t'.join([rv1,rv2,conf,n_genomes])
    out_f.write(line_out + '\n')

out_f.close()
