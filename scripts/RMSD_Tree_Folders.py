#!~/user/bin/env python

print "Creating Tree"

import numpy as np

import pandas as pd

import itertools

from Bio.Phylo.TreeConstruction import DistanceTreeConstructor

from Bio.Phylo.TreeConstruction import _DistanceMatrix

from Bio import Phylo

from Bio import motifs

from Bio.Seq import Seq

from Bio.Alphabet import IUPAC

import copy

import os.path
from os.path import isfile

#Dave modified the matrix.py file so that it includes U in the dictionary

#This creates a function called "walk_the_tree" which takes a path as the input.
def walk_the_tree(path):
    names = [x.name for x in path] #This creates a list of clade names from the path.
    for n in range(len(names)): #This loops through the list of names.
        c=tree.common_ancestor({"name":"%s" %names[n]}) #Each clade name is used to return the corresponding clade object as 'c'
        files=[x.name for x in c.get_terminals()] #A files list is created which contains the terminal clade names from each clade.
        sequences=[s.split('_')[1] for s in files] #A sequences list is created that contains just the sequence string from each filename.
        for sequence in sequences: #Each sequence in the sequence list is made into a sequence object
            Seq(sequence, alphabet=IUPAC.unambiguous_rna)
        m = motifs.create(sequences, alphabet=IUPAC.unambiguous_rna) #A motif object is created from the sequence list
        f.write("%s: %s %s %7d\n" % #Prints the clade name: consensus sequence, degenerate sequence, and # of sequences in the clade
         (c,m.consensus,m.degenerate_consensus,len(sequences)))
        f.write("%s" %m.counts.normalize()) #Prints the normalized count of each nucleotide at each position

# Creates a function called "loop_the_tree" that "walks_the_tree" from each terminal node.
def loop_the_tree(terminals):
    for name in terminals:
        path = tree.trace({"name": "%s" % name}, tree.root)
        f.write("This is the beginning of %s path.\n" % name)
        walk_the_tree(path)

# This creates a function called rename_tree2. The function input is a path.
# The clades in tree2 along that path are renamed as the degenerate sequence for that clade.
def rename_tree2(path):
    names = [x.name for x in path]
    for n in range(len(names)):
        c = tree.common_ancestor({"name": "%s" % names[n]})
        files = [x.name for x in c.get_terminals()]
        sequences = [s.split('_')[1] for s in files]
        for sequence in sequences:
            Seq(sequence, alphabet=IUPAC.unambiguous_rna)
        m = motifs.create(sequences, alphabet=IUPAC.unambiguous_rna)
        newname = tree2.find_any({"name": "%s" % c.name})
        if newname is not None:
            newname.name = "%s" % m.degenerate_consensus

# This creates a function called rename_the_tree. It walks through all of the paths in the original tree and renames tree2.
def rename_the_tree(terminals):
    for name in terminals:
        path = tree.trace({"name": "%s" % name}, tree.root)
        rename_tree2(path)

def check_folder(folder):
    if os.path.isdir(folder) == False:
        print('Folder\'%s\' does not exist...creating' % folder)
        try:
            os.makedirs(folder)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

#Change the 'string' to the filename of interest
filename = 'RMSD.txt'

#Creates a dataframe from the file, the information is tab seperated
df = pd.io.parsers.read_csv(filepath_or_buffer=filename,sep='\t')

# Group the dataframe by Seq1 values
rgroup = df.groupby(df.Seq1.values)

# Creates an array for each group....array ends up being the last group
for r in rgroup:
    array = r[1]

#Generates data from text file, skipping the first line
#Data is a list of lists
data = np.genfromtxt(filename,skip_header=1)

#Delete the nan from the data so that it only contains the RMSD values
data=data[:,2:]

#This flattens the list of lists so that all RMSD values are contained in a single list
flatlist=list(itertools.chain.from_iterable(data))

#This creates a new list that is a list of lists, where each list is the length of the seq1 values in each group
i=0
n=len(array.Seq1.values)
new_list=[]
while i<len(flatlist):
  new_list.append(flatlist[i:i+n])
  i+=n

#This creates a new list of lists called matrix. Each list in the new list is shortened so that it stops after the zero
v=1
matrix=[]
for l in new_list:
    l=l[:v]
    v=v+1
    matrix.append(l)

#This creates an array of all of the Seq2 names in the last group
names_array=array.Seq2.values

#This turns the names array into a list
names=names_array.tolist()

constructor = DistanceTreeConstructor()

#This creates a distance matrix from the names and matrix defined above
dm=_DistanceMatrix(names,matrix)

#Construct tree using unweighted pair group method analysis (UPGMA)
tree = constructor.upgma(dm)

#Sorts clades according to number of terminal nodes
tree.ladderize()

#Write tree to phyloxml file
Phylo.write([tree], 'UPGMAtree.xml', 'phyloxml')

#Read the phyloxml file
tree = Phylo.read('UPGMAtree.xml', 'phyloxml')

#Can delete this line...only use to visualize tree in terminal
#Phylo.draw_ascii(tree)

#Can delete this line...only use to visualize the tree as a string (includes branch lengths)
#print(tree)

#Creates a text file called 30tetraloop_seq_analysis
f=open('seq_analysis.txt', 'w')

#This creates a list of all of the terminal nodes in the tree
terminals=[x.name for x in tree.get_terminals()]

#Execute the function "loop_the_tree" on the terminals list.
loop_the_tree(terminals)

#This creates a complete copy of the tree called tree2, so that tree2 can be modified, while keeping the original tree intact.
tree2 = copy.deepcopy(tree)

#This executes the rename_the_tree function using the terminals list.
rename_the_tree(terminals)

#Write tree2 to phyloxml file
Phylo.write([tree2], 'UPGMAtree_degenerate.xml', 'phyloxml')

#Assign distance cutoff---Half of the RMSD cutoff
RMSD_cutoff = input('RMSD cutoff: ')
distance_cutoff = RMSD_cutoff/2

print "Creating folders"

#Determine the highest-order clade number
count_clades = tree.count_terminals() - 1

#Create a list, sorted ascending, of every clade name
non_terminals = []
i = 1
while i <= count_clades:
    name = "Inner" + str(i)
    non_terminals.append(name)
    i = i + 1

clades_to_collapse= [] #creates list
for clade in non_terminals: #loops through each clade in non-terminals
    clade_name = tree.find_clades({"name": clade}) #returns iterable of clade object with specified name
    terminals_to_collapse = [] #creates list
    for n in clade_name: #for clade object in interable
        terminals=[x.name for x in n.get_terminals()] #creates a list containing the name of all of the terminals for each clade
        for name in terminals: #Loops through the names in the terminals list
            distance=tree.distance({'name':'%s'%name},{'name':'%s'%n}) #Calculates the distance from the terminal clade to the inner clade
            if distance < distance_cutoff: #If distance is less than distance_cutoff
                terminals_to_collapse.append(name) #Add terminal to list
            if set(terminals_to_collapse) == set(terminals): #If all of the terminals in these lists match
                clades_to_collapse.append(n) #Add the clade to the list


#Collapse the clades in the list, preserving the branch length
for clade in clades_to_collapse:
    clade.collapse_all()

#Creates a new iterable that contains the clades that aren't terminal after the collapse
new_non_terminals=tree.find_clades(terminal=False)


for clade in new_non_terminals: #Loop through each clade in new_non_terminals
    terminals=[x.name for x in clade.get_terminals()] #Creates a list containing the names of all of the terminals in each clade
    sequences=[s.split('_')[1] for s in terminals]
    for sequence in sequences:
        Seq(sequence, alphabet=IUPAC.unambiguous_rna)
        m = motifs.create(sequences, alphabet=IUPAC.unambiguous_rna)
        d = m.degenerate_consensus
    terminals_to_folder = [] #Create empty terminals_to_folder list
    for name in terminals: #Loops through each name in terminals list
        distance=tree.distance({'name':'%s'%name},{'name':'%s'%clade}) #Calculates the distance between the terminal clade and inner clade
        if distance < distance_cutoff: #If distance is less than the distance_cutoff
            terminals_to_folder.append(name) #Add terminal clade to terminals_to_folder list
        if set(terminals_to_folder) == set(terminals): #If the terminals_to_folder list equals the terminals list
            folder="%s"%d #folder is equal to the name of the clade
            check_folder(folder) #Run check folder function on folder...see if a folder exists with that clade name
            for name in terminals_to_folder:
                if isfile("%s" %name): #If the name of the terminal is the name of a file
                    start="%s" %name #Starting place is equal to the name of the file
                    destination="%s\%s" %(folder,name) #End destination is equal to the folder/name
                    os.rename(start,destination) #Move the file to the destination folder
