
# coding: utf-8


#!/usr/bin/env python


from os import path
import glob #KR-glob is a module that finds all the pathnames matching a specified pattern according to the rules of Unix shell,
#KR-although the results are returned in arbitrary order
import numpy as np


sequence = '*.pdb' #Set sequence identity



seq_list = glob.glob('%s'%sequence) #Creates a list that contains all pdb files in the path



#Creates a function that removes .pdb from filename and returns the basename
def get_basename(filename):
    
    basename = path.basename(path.splitext(filename)[0])
    return basename


class PDB: #Defines a class that deals with PDB objects

    def __init__(self, filename):
        self.filename   = filename #Returns the filename of a PDB object
        self.basename   = get_basename(filename) #Runs the get_basename function on PDB object (returns filename w/o .pdb)
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

def get_phosphate_coords(fname):
    #returns list of x,y,z, coords of all phosphate atoms in fname
    with open(fname) as f:
        content = f.readlines()
    coords = []
    for line in content:
        record = line[0:6].strip()
        atomname = line[12:16].strip()
        if (record == 'ATOM' or record == 'HETATM' or record == 'HETAT') and "P" in atomname:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coords.append([x,y,z])
    return np.asarray(coords)

def get_backbone_and_sugar_coords(fname):
    #returns list of x,y,z, coords of all backbone and sugar atoms in fname
    with open(fname) as f:
        content = f.readlines()
    coords = []
    for line in content:
        record = line[0:6].strip()
        atomname = line[12:16].strip()
        if (record == 'ATOM' or record == 'HETATM' or record == 'HETAT') and (("P" in atomname) or ("'" in atomname)):
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coords.append([x,y,z])
    return np.asarray(coords)

def get_backbone_sugar_and_selectbase_coords(fname):
    #returns list of x,y,z, coords of all backbone, sugar, and select base atoms in fname
    with open(fname) as f:
        content = f.readlines()
    coords = []
    for line in content:
        record = line[0:6].strip()
        atomname = line[12:16].strip()
        basename = line[19:20]
        if (record == 'ATOM' or record == 'HETATM' or record == 'HETAT') and (((atomname == 'N9' or atomname == 'C8' or atomname == 'C4') and
            (basename == 'G' or basename =='A')) or ((atomname == 'N1' or atomname == 'C2' or atomname == 'C6') and
            (basename == 'C' or basename == 'U')) or ("P" in atomname) or ("'" in atomname)):
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            coords.append([x,y,z])
    return np.asarray(coords)

#Assigns the fixed atoms
#fit_atoms = ["P","O5'","C2'"]
#fit_atoms = ["phosphate"]
#fit_atoms = ["backbone"]
fit_atoms = ["backbone and base"]


results = [] #Creates a list called results

fix_list=seq_list #Sets fix_list as seq_list
mov_list=seq_list #Sets mov_list as seq_list
unique_pdbs = list(set(fix_list).union(set(mov_list))) #Merge the two lists for reading (avoids reading duplicates if in both arrays)

# Read the PDB data and return a dictionary PDB objects
pdbs = read_pdbs(unique_pdbs, fit_atoms)#KR-Dictionary contains the filenames as keys and the pdb objects as values

# Do the calculations
for fix in fix_list: #KR-For each file in fix_list, do the following:
    
    f = pdbs[fix] #KR-Selects the fixed pdb object from pdbs
    fname = f.filename #KR-for each fixed pdb, return the filename
    if fit_atoms[0] == "phosphate":
        coords1 = get_phosphate_coords(fname)
        COM1 = calculate_COM(coords1)
    if fit_atoms[0] == "backbone":
        coords1 = get_backbone_and_sugar_coords(fname)
        COM1 = calculate_COM(coords1)
    if fit_atoms[0] == "backbone and base":
        coords1 = get_backbone_sugar_and_selectbase_coords(fname)
        COM1 = calculate_COM(coords1)

    for mov in mov_list: #KR-For each file in mov_list do the following:

        m = pdbs[mov] #KR-Selects the move pdb object from pdbs
        mname  = m.filename #KR-for each moved file, return the filename
        if fit_atoms[0] == "phosphate":
            coords2 = get_phosphate_coords(mname)
            COM2 = calculate_COM(coords2)
        if fit_atoms[0] == "backbone":
            coords2 = get_backbone_and_sugar_coords(mname)
            COM2 = calculate_COM(coords2)
        if fit_atoms[0] == "backbone and base":
            coords2 = get_backbone_sugar_and_selectbase_coords(mname)
            COM2 = calculate_COM(coords2)

        U, RMSD = calculate_rotation_rmsd(coords1,coords2,COM1,COM2)

        if (U is None) or (RMSD is None):
            continue

        results.append([fname,mname,U,RMSD]) #Adds the fix_filename, mov_filename, U, and RMSD to results list


text_file = open('rmsd.txt',"w") #Creates a text file called rmsd.txt



text_file.write("Seq1\tSeq2\tRMSD\n") #Writes the header "Seq1..tab..Seq2..tab..RMSD" to the text file
for r in results: #Loops through the results list and for each result:
    text_file.write("{:s}\t{:s}\t{:8.3f}\n".format(r[0],r[1],r[3])) #Writes the results to the text file
    #'fix_filename..tab..mov_filename..RMSD' is written to the textfile


#Closes the text file
text_file.close()


