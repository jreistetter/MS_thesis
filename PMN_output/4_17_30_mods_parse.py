#Script to parse the output of the 4.17.12.30_mods output

#Written by Joe Reistetter

class Mod:
    """Holds a module name, its parents, and CPD"""
    def __init__(self, tup):
        self.name = tup[0]
        self.parents = tup[1]
        self.cpd = tup[2]

    def print_mod(self, file_obj):
        if self.parents[0] == '':
            line = self.name + '\t' + '\t'.join(['NA']*5)
            file_obj.write(line + '\n')

        else:
            for cpd_list in self.cpd:
                line = self.name.strip() + '\t'
                line += ' '.join(self.parents)
                line += '\t' + cpd_list[0] + '\t'
                line += '\t'.join(cpd_list[1])
                file_obj.write(line + '\n')

class GeneList:
    """Holds the genes assigned to a module"""
    def __init__(self, tup):
        self.name = tup[0]
        self.genes = tup[1]

    def print_mod(self, file_obj):
        for gene in self.genes:
            line = self.name + '\t' + gene + '\n'
            file_obj.write(line)

class Pathway:
    """Holds the genes that are pathway leading to module"""
    def __init__(self, tup):
        self.name = tup[0]
        self.genes = tup[1]

    def print_mod(self, file_obj):
        line = self.name + ':' + '\t'.join(self.genes) + '\n'
        file_obj.write(line)
        
    

def parse_par_block(par_block):
    vals = par_block.split("'")
    mod_name = vals[1]
    
    #mod_parents looks like '(RV3244C RV0237) '
    mod_parents = vals[2]
    #remove surrounding whitespace
    mod_parents = mod_parents.strip()
    #remove parens and separate the values
    mod_parents = mod_parents[1:-1].split(' ')
    
    #pull out CPD for each combo of parent values
    cpds_split = vals[3].split(') (')
    cpds = []
    for cpd in cpds_split:
        cpd_vals = cpd.split(')')
        cpd_vals[0] = cpd_vals[0].strip('()')
        cpd_vals[1] = cpd_vals[1].split()
        cpds.append(cpd_vals)

    #Last three )) screw up the split, remove the empty list elements
    cpds[-1] = cpds[-1][:2]
    return tuple([mod_name, mod_parents, cpds])

def parse_modules(path, out):
    f = open(path, 'r')
    txt = f.read()
    #Pull out the parent assignments
    #Split on "Scheme", use first portion
    par_raw = txt.split('Scheme')[0]
    #Pull out each block of text representing one module
    #Last element is blank so exclude that as well.
    par_blocks = par_raw.split('\n\n')[:-1]

    out_f = open(out, 'w')
    header = ['moduleID', 'parents', 'parent_state', 'mod_down', 'mod_nochange', 'mod_up']
    out_f.write('\t'.join(header) + '\n')

    mods = []
    for block in par_blocks:
        mod = Mod(parse_par_block(block))
        mod.print_mod(out_f)

def parse_members(path, out):
    f = open(path, 'r')
    txt = f.read()
    mod_txt = txt.split('Assignments.............\n')[1].strip().split('\n')

    out_f = open(out, 'w')
    headers = ['moduleID', 'gene']
    out_f.write('\t'.join(headers) + '\n')
    
    
    for mod_chunk in mod_txt:
        vals = mod_chunk.split(':')
        mod_id = vals[0].split('{')[1].strip(' }')

        mod_members = []

        for val in vals[1].split():
            if 'RV' in val:
                mod_members.append(val)

        mod = GeneList(tuple([mod_id, mod_members]))
        mod.print_mod(out_f)

def parse_pathways(path, out):
    f = open(path, 'r')
    txt = f.read()
    path_chunks = txt.strip().split('\n')
    output = []

    for chunk in path_chunks:
        vals = chunk.split(':')
        mod_id = vals[1].strip().split(' ')[1]

        proteins = []
        
        for val in vals[0].split(' '):
            if 'RV' in val:
                #List has pRV0667, remove p and append
                proteins.append(val[1:])

        output.append(tuple([mod_id, proteins]))

    return output
            
            
   
    
def main():
    import os
    root_path = os.path.expanduser('~/Dropbox/thesis_work/PMN_output/')
    mods_path = root_path + '4.17.12.30_mods_output.txt'
    protein_path = root_path + '4.17.12.30_mods_output_pathways.txt'
    
    out_path = root_path + '4.17.30_mods_parsed.txt'
    members_out = root_path + '4.17.30_mods_members.txt'
    
    
    parse_modules(mods_path, out_path)
    parse_members(mods_path, members_out)


if __name__ == '__main__':
    main()


    
