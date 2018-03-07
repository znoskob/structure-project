#!/usr/bin/env python


import os
from os import path
import glob
import numpy as np
import shutil


# Creates a function that removes .pdb from filename and returns the basename
def get_basename(filename):
    
    basename = path.basename(path.splitext(filename)[0])
    return basename


class PDB:

    def __init__(self, filename):
        self.filename   = filename
        self.basename   = get_basename(filename)
        self.parsed     = False
        self.contents   = None
        self.atomnames  = None
        self.coords     = None
        self.selidx     = None
        self.COM        = None
        
        self.read_file()
        
        return None
        
    def set_selidx(self,select=None):
        
        if (select == None):
            self.selidx = np.arange(0,len(self.coords))
            return
            
        t = np.asarray(select).reshape(-1,1)
        q = self.atomnames
        
        selidx = np.sort(np.where(q == t)[1])
        self.selidx = selidx
        
        return
        
    def set_coords(self,select=None):
        
        if (self.parsed == False):
            self.read_file()
        
        coords = []
        names  = []
        for line in self.contents:
            
            record = line[0:6].strip()
            if record == 'ATOM' or record == 'HETATM' or record == 'HETAT':
            
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            
                coords.append([x,y,z])
                names.append(line[12:17].strip())
        
        self.coords    = np.asarray(coords)
        self.atomnames = np.asarray(names)
        
        if select != None:
            self.set_selidx(select)
        
        return
    
    def read_file(self):
        
        f = open(self.filename,'r')
        self.contents = f.read().splitlines()
        self.parsed = True
        f.close()
        
        self.set_coords()
        
        return
        
    def get_selcoords(self):
        coords = self.coords[self.selidx]
        return calculate_COM(coords), coords



def calculate_COM(coords):

    L = coords.shape[0]
    COM = np.sum(coords,axis=0) / float(L)
    
    return COM


def read_pdbs(filelist,fit_atoms=None):
    
    results = {}
    
    for f in filelist:
        
        pdbf = PDB(f)
        COM = calculate_COM(pdbf.coords)
        pdbf.COM = COM
        
        if fit_atoms != None:
            pdbf.set_selidx(fit_atoms)
            
        results[f] = pdbf 
        
    return results



def calculate_rotation_rmsd(coords1,coords2,COM1,COM2):

    sel1 = coords1 - COM1
    sel2 = coords2 - COM2

    # check for consistency
    if len(sel1) != len(sel2):
        return None, None
        
    L = len(sel1)
    assert L > 0


    # Initial residual, see Kabsch.
    R0 = np.sum( np.sum(sel1 * sel1,axis=0),axis=0) + np.sum( np.sum(sel2 * sel2,axis=0),axis=0)

    # Calculate the components of the rotation matrix (V,W)
    # S is used to calculate the error (RMSD)
    V, S, W = np.linalg.svd(np.dot( sel2.T, sel1))

    # Calculate if the poduct of the determinants is + or -
    # if negative reflect the rotation matrix components prior
    # determining the rotation matrix (U)
    reflect = float(str(float(np.linalg.det(V) * np.linalg.det(W))))

    if reflect == -1.0:
    	S[-1] = -S[-1]
    	V[:,-1] = -V[:,-1]

    U = np.dot(V, W)

    # Calculate the RMSD using sigma from the SVD calculation above
    RMSD = R0 - (2.0 * sum(S))
    RMSD = np.sqrt(abs(RMSD / L))

    return U, RMSD



def translate_rotate_coords(coords,COM,U=None):
    
    # Translate only
    if U == None:
        return coords - COM
        
    # Translate and rotate
    return np.dot((coords-COM),U)

# This is a function to remove redundancies from list
def make_unique(original_list):
    unique_list = []
    [unique_list.append(obj) for obj in original_list if obj not in unique_list]
    return unique_list

def atom_types_in_file(fname,unique):
    #returns list of atoms in order from file when unique=0
    #returns list of unique atoms in file when unuque=1, assume random order
    atoms = []
    with open(fname) as f:
        content = f.readlines()
    for line in content :
        record = line[0:6].strip()
        if record == 'ATOM' or record == "HETATM" or record == 'HETAT':
            atoms.append(line[13:16].strip())
    if unique == 0:
        return atoms
    else:
        return list(set(atoms))    
  
def atom_count_in_file(fname):
    #returns the number of atoms in the file
    with open(fname) as f:
        content = f.readlines()        
    return (int(content[len(content)-3][7:11]))  
    
def get_all_coords(fname):
    #returns list of x,y,z, coords of all atoms in fname
    with open(fname) as f:
        content = f.readlines()
    coords = []  
    for line in content:
        record = line[0:6].strip()
        if record == 'ATOM' or record == "HETATM" or record == 'HETAT':
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip()) 
            coords.append([x,y,z])
    return np.asarray(coords)    

#Defines function that gets the backbone and select base coordinates
def get_backbone_sugar_and_selectbase_coords(fname):
    # returns list of x,y,z, coords of all backbone, sugar, and select base atoms in fname
    with open(fname) as f:
        content = f.readlines()
    coords = []
    for line in content:
        record = line[0:6].strip()
        atomname = line[12:16].strip()
        basename = line[19:20]
        if (record == 'ATOM' or record == 'HETATM' or record == 'HETAT') and (
                    ((atomname == 'N9' or atomname == 'C8' or atomname == 'C4') and
                         (basename == 'G' or basename == 'A')) or (
                    (atomname == 'N1' or atomname == 'C2' or atomname == 'C6') and
                    (basename == 'C' or basename == 'U')) or ("P" in atomname) or ("'" in atomname)):
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coords.append([x, y, z])
    return np.asarray(coords)

def atom_types_in_backbone_select_base(fname,unique):
    #returns list of atoms in order from file when unique=0
    #returns list of unique atoms in file when unique=1, assume ransom order
    atoms = []
    print fname
    with open(fname) as f:
        content = f.readlines()
    for line in content :
        record = line[0:6].strip()
        atomname = line[12:16].strip()
        basename = line[19:20]
        if (record == 'ATOM' or record == 'HETATM' or record == 'HETAT') and (
                            ((atomname == 'N9' or atomname == 'C8' or atomname == 'C4') and
                                 (basename == 'G' or basename == 'A')) or (
                                    (atomname == 'N1' or atomname == 'C2' or atomname == 'C6') and
                                    (basename == 'C' or basename == 'U')) or ("P" in atomname) or ("'" in atomname)):
            atoms.append(line[13:16].strip())
    if unique == 0:
        return atoms
    else:
        return list(set(atoms))

def write_all_translated_coords(out_data_path,out_file,in_data_path,template_file,all_avg_coords):
    #creates file seq_avg.pdb and writes all_avg_coords  
    #using template_file
    fp = open(template_file, "r")
    content = fp.readlines()
    fp.close()
    fp = open("new_template.txt", "w+")
    # Only keeps the lines that contain the selected atoms and saves into new template file
    for line in content:
        record = line[0:6].strip()
        atomname = line[12:16].strip()
        basename = line[19:20]
        if (record == 'ATOM' or record == 'HETATM' or record == 'HETAT') and (
                            ((atomname == 'N9' or atomname == 'C8' or atomname == 'C4') and
                                 (basename == 'G' or basename == 'A')) or (
                                    (atomname == 'N1' or atomname == 'C2' or atomname == 'C6') and
                                    (basename == 'C' or basename == 'U')) or ("P" in atomname) or ("'" in atomname)):
            fp.write(line)
    fp.close()
    fp = open("new_template.txt", "r")
    lines = fp.read()
    linesl = lines.splitlines()
      
    avg_file = []
    for i in range(len(all_avg_coords)):
        avg_file.append(linesl[i][0:30]+ \
                  "{:8.3f}".format(all_avg_coords[i][0])+ \
                  "{:8.3f}".format(all_avg_coords[i][1])+ \
                  "{:8.3f}".format(all_avg_coords[i][2])+ \
                  linesl[i][54:])
    avg_file.append("TER")
    avg_file.append("END")
    fp=open(out_data_path+out_file,"w")
    for l in avg_file:
        fp.write(l+"\n")
    fp.close()
    return        


  
def find_rmsd(fname1,fname2,fit_atoms):
 
      #input two files and atoms to fit
      #assumes first file is fixed, second is move atoms
      #returns U,rmsd, and 
      #coords of fit_atoms and all atoms
      #check atom count in both files      
      if atom_count_in_file(fname1) != atom_count_in_file(fname2) :
          print("atom counts do not match: ")
          print(fname1,": ",atom_count_in_file(fname1))
          print(fname2,": ",atom_count_in_file(fname2))
          quit()
      
      #check atom types amd order in both files
      if atom_types_in_file(fname1,0) != atom_types_in_file(fname2,0):
          print('atom types/order do not match: ',fname1," ",fname2)
          quit()
      
      unique_pdbs=[fname1,fname2]
      pdbs = read_pdbs(unique_pdbs,fit_atoms); 
      
      f1 = pdbs[fname1]
      fname1 = f1.basename #KR-for each fixed pdb, return the basename
      COM1, coords1 = f1.get_selcoords()
           
      f2 = pdbs[fname2]
      fname2  = f2.basename #KR-for each moved file, return the basename
      COM2, coords2 = f2.get_selcoords() 
 
      U, rmsd = calculate_rotation_rmsd(coords1,coords2,COM1,COM2)
 
      return U,rmsd 
      
def mark_min_rmsd(results,seq): 
    # enters 1 in col=0 of results for minimum
    #value of rmsd(col=3) associated with 
    #file name (col=2).  Assumes that current
    #col=0 content is zero
    #find rows with correct fname and make col=0 to be -1
    look = []
    val = []
    for i in range(len(results)):
            look.append(i)
            val.append(results[i][4])

    results[look[val.index(min(val))]][0] = 1
    results[look[val.index(min(val))]][1] = len(look)

    return results

def write_results(out_data_path,results_file,seq):
    #writes results values into file
    fp=open(out_data_path+results_file,"w")
    for r in results:
        out_string="{0:3d},{1:3d},  {2:s},{3:7.4f},{4:7.4f}".format(r[0],r[1],r[3],r[4],r[5])
        fp.write(out_string+"\n")
    fp.close()
    return

def purge_bad_files(unique_list):
    #cleans unique_list and removes them from the list when "bad"/ 
    #current definition of bad is (1) when file has only one line
    for un in unique_list:
        for fname in un:
            f = open(fname,'r')
            contents = f.read().splitlines()
            f.close()
            if len(contents) <= 1:
                un.remove(fname)  
                print(fname+" removed")
        if len(un) == 0:        
            unique_list.remove(un)
    return unique_list            
            
    
def purge_h_atoms(unique_list):
    #for all files in unique_list the following is done:
    #(1) copy of file is made in same directory as *-H.pdb
    #(2) H's are removed from existing file and atom num bers are updated.
    for un in unique_list:
        for fname in un:
            atoms = atom_types_in_file(fname,1)
            h_yes = 0
            if "H" in str(atoms):
                h_yes = 1
                print('H in '+fname)
                #copy file
                shutil.copy(fname,fname.split(".pdb")[0]+"_H.pdb")
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
              
    return        
        
def reorder_list(fnames):
    #orders list of file names to have the smallest number of
    #atoms as first
    size=[]
    for l in fnames:
        size.append(atom_count_in_file(l))
    small_atoms_loc = size.index(min(size))
    temp = fnames[0]
    fnames[0] = fnames[small_atoms_loc]
    fnames[small_atoms_loc] = temp
    
    return fnames
    
    
        
#set data path relative to dir of script            
in_data_path = ""
out_data_path = in_data_path

sequence = '*.pdb'

seq_list = glob.glob(in_data_path+'%s'%sequence) # Creates a list that contains all pdb files in the path

directory = os.getcwd()


degseq = os.path.basename(directory)

fit_atoms = ["backbone and select base"]

#use only first file for fix (=1) or cycle through them all (=0)
fix_only_first = 1

results = [] # Creates a blank list called results
rep_structure = [] # This creates a blank list called rep-structure

fix_list = list(seq_list[:1])
mov_list = list(seq_list[1:])
unique_pdbs = list(set(fix_list).union(set(mov_list)))

all_avg_coords = []
avg_coords = []

translated_coords = [0.0,0.0,0.0]
all_translated_coords = [0.0, 0.0, 0.0]
#open file for rmsd output
rmsd_fp = open(out_data_path+degseq+"_rmsd.txt","w")
        
# Reads the pdbs from the unique_pdbs list and creates pdb objects that contain the fit atoms
pdbs = read_pdbs(unique_pdbs,fit_atoms)
        
#included tracks files which are included in calculation
included = []


for fix in fix_list:
    sel1 = []
    f = pdbs[fix]
    fname = f.filename
    fseq = fname.split('_')[1]
    print fit_atoms[0]
    if fit_atoms[0] == "backbone and select base":
        coords1 = get_backbone_sugar_and_selectbase_coords(fname)
        COM1 = calculate_COM(coords1)
    if fit_atoms[0] == "all":
        COM1, coords1 = f.get_selcoords()
    sel1 = coords1 - COM1
    # Creates a blank array with the rows equal to the length of coords and 3 columns
    translated_coords = translated_coords + (np.array(sel1))
    #get all atoms in this file and do transformation
    all_translated_coords = all_translated_coords + (coords1-COM1)
    included.append(1)
            #open file for rmsd value output and enter fixed

    out_string="\n(fix){0:s},{1:7.4f}\n".format(fname,0.0)
    rmsd_fp.write(out_string)

    #mov_list.remove(fix)
    for mov in mov_list:
        sel2 = []
        m = pdbs[mov]
        mname  = m.filename
        if fit_atoms[0] == "backbone and select base":
            coords2 = get_backbone_sugar_and_selectbase_coords(mname)
            COM2 = calculate_COM(coords2)
        if fit_atoms[0] == "all":
            COM2, coords2 = m.get_selcoords()
    
        #coords only added to average calculation if above
        #conditions are met                    
        U, RMSD = calculate_rotation_rmsd(coords1,coords2,COM1,COM2)
        out_string="{0:s},{1:7.4f}\n".format(mname,RMSD)
        rmsd_fp.write(out_string)
        # This uses the rotation matrix U and the COM2 to translate,
        # then rotate the coordinates of the moved structure
        translated_coords2 = translate_rotate_coords(coords2,COM2,U)
        translated_coords = np.array(translated_coords) + np.array(translated_coords2)
        #perform operations on all
        if fit_atoms[0] == "all":
            all_coords=get_all_coords(mname)
                    #translated_coords = translate_rotate_coords(all_coords,COM2,U)
        if fit_atoms[0] == "backbone and select base":
            all_coords=get_backbone_sugar_and_selectbase_coords(mname)
        all_translated_coords = np.array(all_translated_coords) + np.array(translate_rotate_coords(all_coords,COM2,U))
        included.append(1)
                    
rmsd_fp.close()
# This finds the average coordinates by adding the fixed coords to the moved coords and dividing by the
# length of the mov_list plus one for the fix_list
avg_coords = (np.array(translated_coords))/float(sum(included))
all_avg_coords = (np.array(all_translated_coords))/float(sum(included))
        
#write all_avg_coords to file in pdb format
all_avg_fname=degseq+"_avg.pdb"
write_all_translated_coords(out_data_path,all_avg_fname,in_data_path,fname,all_avg_coords)

if fit_atoms[0] == "backbone and select base":
    coords1 = get_backbone_sugar_and_selectbase_coords(all_avg_fname)
    COM1 = calculate_COM(coords1)
sel1 = coords1 - COM1

new_mov_list = seq_list # The new move list is equal to all of the structures in the original list

for mov in new_mov_list: # This loops through the new_mov_list
    if included[new_mov_list.index(mov)] == 1:
        n = pdbs[mov] # This selects the pdb object
        nname = n.filename # This selects the pdb filename
        if fit_atoms[0] == "backbone and select base":
            coords2 = get_backbone_sugar_and_selectbase_coords(nname)
            COM2 = calculate_COM(coords2)
        sel2 = coords2 - COM2
        U, RMSD = calculate_rotation_rmsd(coords1, coords2, COM1, COM2)
        results.append([0,0,all_avg_fname,nname,RMSD,0.0])
#denote minimum for this collection of mov files
#with 1 in col=0
results = mark_min_rmsd(results,degseq)


        
#print to screen
results.sort(key=lambda x:x[4])
for r in results:
    if r[0] == 1:
        print(r)

final_fix_list=[]
final_mov_list=[]
for r in results:
    if r[0] == 1:
        final_fix_list.append(r[3])
    if r[0] == 0:
        final_mov_list.append(r[3])

for fix in final_fix_list:
    sel1 = []
    f = pdbs[fix]
    fname = f.filename
    fseq = fname.split('_')[1]
    if fit_atoms[0] == "backbone and select base":
        coords1 = get_backbone_sugar_and_selectbase_coords(fname)
        COM1 = calculate_COM(coords1)
    if fit_atoms[0] == "all":
        COM1, coords1 = f.get_selcoords()
    sel1 = coords1 - COM1


    for mov in final_mov_list:
        sel2 = []
        m = pdbs[mov]
        mname = m.filename
        if fit_atoms[0] == "backbone and select base":
            coords2 = get_backbone_sugar_and_selectbase_coords(mname)
            COM2 = calculate_COM(coords2)
        if fit_atoms[0] == "all":
            COM2, coords2 = m.get_selcoords()


        U, RMSD = calculate_rotation_rmsd(coords1, coords2, COM1, COM2)

        for r in results:
            if r[3] == mov:
                r[5] = RMSD


#write detailed retsults to file
results_file=degseq+"_results.txt"
write_results(out_data_path,results_file,degseq)

#copy best structures into results folder
if os.path.isdir(out_data_path+"representative"):
    os.rename(out_data_path+"representative",out_data_path+"representative2")
    shutil.rmtree(out_data_path+"representative2")
os.makedirs(out_data_path+"representative")
for r in results:
    if r[0] == 1:
        shutil.copy(r[3],out_data_path+"representative\\")
                

        