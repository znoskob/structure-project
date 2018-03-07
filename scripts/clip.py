#!/usr/bin/env python

from subprocess import call

import shutil
import os.path
from os.path import isfile

from Bio.PDB import *
from optparse import OptionParser

class RangeSelect(Select):
    def __init__(self, chain, resids2, positions_resids2_dict):
        self.chain = chain
        self.resids2 = resids2
    def accept_residue(self, residue):
        if residue.parent.id == self.chain:
            rid = residue.get_id()
            rid = str(rid[1]) + rid[2]
            rid = rid.strip()
            found = 0
            if rid in resids2:
                for key, value in positions_resids2_dict.iteritems():
                    if rid == value:
                        residue.id = (rid[0], key, ' ')
                found = 1
            return found
            
def read_csv(file):
    csvfile = open(file,'r')
    
    #replace csv library processing
    #if list of files was opened with excel it is 
    #possible that some lines start/stop with "
    #and then contents are not split correctly.
    lines = [line.strip() for line in csvfile]
    csvfile.close()
    
    header = lines[0].strip('"').split('#')
    data = []
    for l in lines[1:]:
        data.append(l.strip('"').split('#'))
      
    return header, data

def check_folder(folder):
    if os.path.isdir(folder) == False:
        print('Folder \'%s\' does not exist...creating' % folder)
        try:
            os.makedirs(folder)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

def in_order(resids, beg, end):
    #determines if resid velus are in increments of one or not
    if (end - beg + 1) == (len(resids.split())):
        return 1
    else:
        return 0

def write_rejects(fname,list_name):
    f = open(fname,"w")
    for line in list_name:
        f.write(line+"\n")
    f.close()

def atom_types_in_file(fname,unique):
    #returns list of atoms in order from file when unique=0
    #returns list of unique atoms in file when unuque=1, assume random order
    atoms = []
    with open(fname) as f:
        content = f.readlines()
    for line in content :
        record = line[0:6].strip()
        if record == 'ATOM' or record == "HETATM":
            atoms.append(line[13:16].strip())
            #print(line[13:16])        
    if unique == 0:
        return atoms
    else:
        return list(set(atoms)) 

def purge_h_atoms(fname,h_removed_dir):
    #if H-atoms in file
    #(1) copy of file is made in same directory as *-H.pdb
    #(2) H's are removed from existing file and atom numbers are updated.
    base_fname = os.path.basename(fname)
    atoms = atom_types_in_file(fname,1)
    h_yes = 0
    if "H" in str(atoms):
            h_yes = 1
            print('H in '+fname)
            #copy file
            shutil.copy(fname,h_removed_dir+base_fname.split(".pdb")[0]+"_H.pdb")
            #read file with H
            f = open(fname,'r')
            old_contents = f.read().splitlines()
            f.close()
            anumber=1
            new_contents = []
            for ln in old_contents:
                if "H" in ln[12:16] :
                    pass
                else:
                    out_string="{}{:6d}{}".format(ln[0:5],anumber,ln[11:])
                    new_contents.append(out_string)
                    anumber = anumber +1
    if h_yes == 1:
        f = open(fname,"w")
        for ln in new_contents:
            f.write(ln+"\n")
        f.close()
    return h_yes

# Parse the options and setup input parameters
optparser = OptionParser()

optparser.add_option('--prepend', default=None, dest='prepend', help='number of closing bases to prepend to the sequence [default: 0]')
optparser.add_option('--append', default=None, dest='append', help='number of closing bases to append to the sequence [default: 0]')
optparser.add_option('--overwrite', action="store_true", default=False, dest='overwrite', help='Overwrite existing clipped files if specified, otherwise only missing ones will be created')

(options, args) = optparser.parse_args()


append = 0
prepend = 0
overwrite = options.overwrite
#assign file name or read from user
csvfile = args[0]
#csvfile = "Tetraloops.csv"

if options.prepend is not None:
    print('Prepending %s closing bases' % options.prepend)
    prepend = int(options.prepend)
if options.append is not None:
    print('Appending %s closing bases' % options.append)
    append = int(options.append)

data_path = ".\\"  #path from program to data file
pdb_ext = ".PDB"
pdb_dir = data_path+"pdbs\\"
clips_dir = data_path+"clips\\"
h_removed_dir = clips_dir+"h_removed\\"


# Housekeeping... make sure output directories exist, if not create them
check_folder(pdb_dir)
check_folder(clips_dir)
check_folder(h_removed_dir)

#minimum number of lines in pdb file to be valid
min_pdb_lines = 20

# Read the CSV file. This uses the default CoSSMos delimiter 
# which is the '#' symbol
data = []
header, ldata = read_csv(csvfile)
data.extend(ldata)

pdbidx = header.index('pdb')
seqidx = header.index('Aseq')
if 'Bseq' in header:
    seqBidx = header.index('Bseq')
asnidx = header.index('Aseq_num')
if 'Bseq_num' in header:
    bsnidx = header.index('Bseq_num')
conidx = header.index('A Residue Conformation 2')
expidx = header.index('experiment')
residx = header.index('resolution')

#download all pdb files.
pdb_names = []
for i in range(len(data)):
    data_row = data[i]
    pdbid = data_row[pdbidx]
    
    print("Downloading pdb id: "+pdbid)
    
    # Check for the full PDB file and download if needed
    # Hack: The file gets downloaded to the current directory and
    #       then moved to the pdbs folder
    pdbfile = "%s.pdb" % pdbid
    pdbname = pdb_dir+pdbfile
    pdb_names.append(pdbname)
    if isfile(pdbname) == False:
        print('PDB file %s does not exist...downloading' % pdbid)
        url =  "http://www.rcsb.org/pdb/files/%s.pdb" % pdbid
        call(["curl","-L","-O","-s",url])
        if isfile(pdbfile):
            os.rename(pdbfile,pdbname)
    else:
        print('PDB file %s already exists' % pdbid)
        
#check all pdb files and be sure there is content.
unique_pdb_names = list(set(pdb_names))        
pdbfile_rejects = []        
for nm in unique_pdb_names:
    count = len(open(nm).readlines(  ))
    if count < min_pdb_lines:
        pdbfile_rejects.append(nm)
    else:
        print('PDB file %s size OK ' % nm)

#begin clipping loop
period_rejects = []
h_atoms_removed = []
parser = PDBParser(QUIET=True)
for i in range(len(data)):
    data_row = data[i]
    pdbid = data_row[pdbidx]
    
    print("Loop Processing pdb id: "+pdbid)
    
    pdbfile = "%s.pdb" % pdbid
    pdbname = 'pdbs/%s' % pdbfile
    if pdbname in pdbfile_rejects:
        continue

    chain_id = data_row[conidx].split()[0].rstrip('1234567890').replace("'","")[0]
    if len(chain_id) == 0:
        chain_id = " "

    first_base = data_row[conidx].split()[0]

    if 'Bseq' in header:
        basename = "%s_%s_%s" % (pdbid,data_row[seqidx] + data_row[seqBidx],first_base)
    else:
        basename = "%s_%s_%s" % (pdbid, data_row[seqidx], first_base)
    clipname = clips_dir+basename+pdb_ext
                
    structure = parser.get_structure('asdf',pdbname)
    model = structure[0]
    chain = model[chain_id]
    if 'Bseq' in header:
        resids = data_row[asnidx] + data_row[bsnidx]
    else:
        resids = data_row[asnidx]
        
    
    if '.' in resids:
        print('Contains period.  Deleting periods in %s' % clipname)
        period_rejects.append(clipname)
        print resids
        resids = resids.replace('.','')
        print resids

    
    print('Processing sequence: %s' % basename)
    

    if resids[0] == '-':
        print ('Negative resid')
        resids_split = resids.split()
        resids2 = []
        posresids = []
        for r in resids_split:
            r2 = r.rsplit('-')
            posresids.append(r2[1])
        for r in posresids:
            r2 = "-" + r
            resids2.append(r2)
        positions_resids2 = zip(range(1, len(resids2) + 1), resids2)
        positions_resids2_dict = dict(positions_resids2)

    else:
        resids_split = resids.split()
        resids2 = []
        for r in resids_split:
            r2 = r.rsplit('-')
            resids2.append(r2[0])
        positions_resids2 = zip(range(1, len(resids2) + 1), resids2)
        positions_resids2_dict = dict(positions_resids2)

    # RangeSelect has been modified to renumber the residues
    # matching the selection criteria, starting at 1
    rselect = RangeSelect(chain_id, resids2, positions_resids2_dict)
    
    # Save out the PDB file with modified residue numbers
    io = PDBIO()
    io.set_structure(structure)
    io.save(clipname,rselect)

        
    removed = purge_h_atoms(clipname,h_removed_dir)
    if removed == 1:
        print ("H-atoms removed from "+clipname)
        h_atoms_removed.append(clipname)
    else:
        print ("no H-atoms in "+clipname)
        


print('All done processing CSV file: %s' % csvfile)
#print changes and info to files to clips dir
write_rejects(clips_dir+"pdbfile_rejects.txt",pdbfile_rejects)
write_rejects(clips_dir+"periods_removed.txt",period_rejects)
write_rejects(clips_dir+"H_atoms_removed.txt",h_atoms_removed)
    








