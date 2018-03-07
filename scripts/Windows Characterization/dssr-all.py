# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 20:28:38 2017

@author: kirkpacc
"""
import json
import glob
import os
import subprocess
import sys

def clean_json (fname):
    f1 = open(fname, 'r')
    f2 = open(fname+".tmp", 'w')
    for line in f1:
        f2.write(line.replace('\\', '\\\\'))
    f1.close()
    f2.close()
    os.remove(fname)    
    os.rename(fname+".tmp",fname)

def find_it (look_here,ind1,fld1,ind2,fld2):
    'finds value'
    if (fld1 in look_here[ind1]) and (fld2 in look_here[ind1][fld1][ind2]) :
        return str(look_here[ind1][fld1][ind2][fld2])
    else:
        return "blank"   
        
def dict_count (dict,elem):
    'returns num of occurances of elem in dict'
    if elem in dict:
        return  range(len(dict[elem])) 
    else:
        return range(0,0);

#determine if version 2 or 3
if sys.version_info < (3,0,0):
    pversion = 2
    import Tkinter, tkFileDialog
    from Tkinter import *
else:
    pversion = 3
    import Tkinter, tkFileDialog
    from tkinter import *    


root = Tk()
root.withdraw()
#root.deiconify()
#root.destroy()
#root.fileName = filedialog.askopenfilename(filetypes=(("text","*.txt"),("All","*.*")))
#print(root.fileName)

#path to DSSR executable
if pversion == 2:
    root.fileName = tkFileDialog.askopenfilename(title="Select x3dna-dssr.exe",filetypes=(("exe","*.exe"),("All","*.*")))
    print(root.fileName)
else:
    root.fileName = tkFileDialog.askopenfilename(title="Select x3dna-dssr.exe",filetypes=(("exe","*.exe"),("All","*.*")))
    print(root.fileName)    
DSSR_run_path = root.fileName

#path to dssr input/output files
#output to same directorry as PDB files
if pversion == 2:
    root.dirName = tkFileDialog.askdirectory(title="Choose directory of PDB files")
else:
    root.dirName = tkFileDialog.askdirectory(title="Choose directory of PDB files")        
#print root.dirName
#DSSR_data_path = ".\\"
DSSR_data_path = root.dirName+"\\"
#print(DSSR_data_path)
all_files = glob.glob(DSSR_data_path+"*")
pdb_files = []
for f in all_files:
    if 'PDB' in f:
        pdb_files.append(f);

#run DSSR for each PDB file
for file in pdb_files:
    basefname = file.replace(".PDB","")
    args = " --input="+file+" --json -o="+basefname+".json --prefix=temp --more --non-pair --u-turn --po4 --idstr=long "
    run_string = DSSR_run_path+args
    subprocess.call(run_string)
    #json file needs to be cleaned of single escape characters
    clean_json (basefname+".json")

#delete merged.json if it exists
if os.path.isfile(DSSR_data_path+"merged.json"):
    os.remove(DSSR_data_path+"merged.json")
    #print("File Removed!")

#list files
jsons = glob.glob(DSSR_data_path+'*.json')
#print (jsons)

#loop through files--append contents and ,
f = open(DSSR_data_path+"merged.json", "w")
f.write("[\n")
for fname in jsons:
    fr = open(fname)
    f.write(fr.read())
    if (fname != jsons[-1]): 
        f.write(",\n")
f.write("]")
f.close()
        
#delete all files with form temp-*
for f in glob.glob(DSSR_data_path+"temp-*"):
    os.remove(f)    

# code for json_bp.py


#from tkinter import filedialog
#from tkinter import *

#output file name
fname=".\\bp.txt"
#following outputs all data to file when =1, 
#smaller output when =0
out_more=0 

#output to screen (=1) or silent (0)
out_screen=1

#field delimiter and eol
delim="#"
eol="\n"
 
#open json file  
json_fname = ".\\merged.json"
clean_json(json_fname)       
with open(json_fname) as json_file:
    json_data = json.load(json_file) 


#names for column labels. should end with eol character
#User must be sure this list is consistent with
#value of out_more
labels = ['sequence','nt1','nt2','bp','name'] 

#open output file
fp = open(fname, "w")
#write dolumn labels
for l in labels:
    fp.write(l+delim)
fp.write(eol)
    

results = []
for ind in range(len(json_data)):
        f1=json_data[ind]['metadata']['str_id']+delim
        
        for ind2 in dict_count(json_data[ind],'pairs'):
            f2=find_it(json_data,ind,'pairs',ind2,'nt1')
            f2s=f2[len(f2)-2:len(f2)-1]+delim 
            f2l=f2+delim+f2s
              
            f3=find_it(json_data,ind,'pairs',ind2,'nt2')
            f3s=f3[len(f3)-2:len(f3)-1]+delim
            f3l=f3+delim+f3s
                            
                            
            f4=find_it(json_data,ind,'pairs',ind2,'bp')
            f4s=f4[len(f4)-1:len(f4)-1]+delim
            f4l=f4+f4s               
                            
            f5=find_it(json_data,ind,'pairs',ind2,'name')             
        
            

            if out_more == 1:
                fp.write(f1+f2l+f3l+f4l+f5+eol)
                if out_screen == 1:
                    print(f1+f2l+f3l+f4l+f5)
            else:
                fp.write(f1+f2s+f3s+f4l+f5+eol)
                if out_screen == 1:
                    print(f1+f2s+f3s+f4l+f5)

fp.close()
   
#code for json_stacks.py
   
def find_in_str(str):
    #returns the single occurence of a string in list that is in str
    #exits if none of more than one string in list is found in str 
    #list is given below
    to_look_for=['<<','<>','><','>>']
    x=[int(str.find(to_look_for[0])),int(str.find(to_look_for[1])),int(str.find(to_look_for[2])),int(str.find(to_look_for[3]))]
    if max((x))-3 == sum((x)) :
        pos=x.index(max(x)); 
        return to_look_for[pos]
    else:
        print ('not ok in find_in_str')
        exit()

#output file name
fname=".\\stacks.txt"
#following outputs all data to file when =1, 
#smaller output when =0
out_more=0 

#output to screen (=1) or silent (0)
out_screen=1

#field delimiter and eol
delim="#"
eol="\n"


#json_fname = '.\\merged.json'
json_fname = ".\\merged.json"
clean_json(json_fname)  
with open(json_fname) as json_file:
    json_data = json.load(json_file)   

#names for column labels. should end with eol character
#User must be sure this list is consistent with
#value of out_more
labels = ['sequence','nt1','nt2','stacking'] 

 
#open output file
fp = open(fname, "w")
#write dolumn labels
for l in labels:
    fp.write(l+delim)
fp.write(eol)
 
#begin program   
for ind in range(len(json_data)):
    #find sequence
    f1=json_data[ind]['metadata']['str_id']+delim
    
    #determine number of nonPairs. Loop to find nt1,nt2,stacking
    for ind2 in dict_count(json_data[ind],'nonPairs'):
        
        #f2 is nt1
        f2=find_it(json_data,ind,'nonPairs',ind2,'nt1')
        f2s = f2[len(f2)-2:len(f2)-1]+delim   
        f2l=f2+f2[len(f2)-2:len(f2)-1]+delim
        
        #f3 is nt2
        f3=find_it(json_data,ind,'nonPairs',ind2,'nt2')
        f3s = f3[len(f3)-2:len(f3)-1]+delim 
        f3l=f3+f3[len(f3)-2:len(f3)-1]+delim
        
        #f4 is stacking
        f4=find_it(json_data,ind,'nonPairs',ind2,'stacking')+delim
        f4s = find_in_str(f4) 
        f4l = f4+ find_in_str(f4) +delim
        
        outstringl=f1+f2+f3+f4
        
        #output to file only if all fields are not blank
        
        if out_more == 1:
            fp.write(f1+f2l+f3l+f4l+eol)
            if out_screen == 1:
                print(f1+f2l+f3l+f4l)
        else:
            fp.write(f1+f2s+f3s+f4s+eol)
            if out_screen == 1:
                print(f1+f2s+f3s+f4s)
                    
fp.close()  

#code for json_hbond

#output file name
fname=".\\hbonds.txt"
#following outputs all data to file when =1, 
#smaller output when =0
out_more=0 

#output to screen (=1) or silent (0)
out_screen=1

#field delimiter and eol
delim="#"
eol="\n"
 
#open json file  
json_fname = ".\\merged.json"
clean_json(json_fname)       
with open(json_fname) as json_file:
    json_data = json.load(json_file) 

#names for column labels. should end with eol character
#User must be sure this list is consistent with
#value of out_more
labels = ['sequence','nt1','nt2','hbond1','hbond2','hbond3', 'hbond4'] 

#open output file
fp = open(fname, "w")
#write dolumn labels
for l in labels:
    fp.write(l+delim)
fp.write(eol)
    

results = []
for ind in range(len(json_data)):
        f1=json_data[ind]['metadata']['str_id']+delim
        
        for ind2 in dict_count(json_data[ind],'nonPairs'):
            f2=find_it(json_data,ind,'nonPairs',ind2,'nt1')
            f2s=f2[len(f2)-2:len(f2)-1]+delim 
            f2l=f2+delim+f2s
              
            f3=find_it(json_data,ind,'nonPairs',ind2,'nt2')
            f3s=f3[len(f3)-2:len(f3)-1]+delim
            f3l=f3+delim+f3s
                            
            f4=find_it(json_data,ind,'nonPairs',ind2,'hbonds_desc')+delim               
            #remove [  ] 
            f4=re.sub(r'\[.*?\]','', f4)
            #break at comma into two strings
            f4=re.sub(',', delim, f4)
            
            if out_more == 1:
                fp.write(f1+f2l+f3l+f4+eol)
                if out_screen == 1:
                    print(f1+f2l+f3l+f4)
            
            else:
                fp.write(f1+f2s+f3s+f4+eol)
                if out_screen == 1:
                    print(f1+f2s+f3s+f4)



fp.close() 

#code for json_conf_pucker.py

#output file name
fname=".\\conf_pucker.txt"
#following outputs all data to file when =1, 
#smaller output when =0
out_more=0 

#output to screen (=1) or silent (0)
out_screen=1

#field delimiter and eol
delim="#"
eol="\n"

       
#open json file     
json_fname = ".\\merged.json"    
with open(json_fname) as json_file:
    json_data = json.load(json_file)  
    

#names for column labels. should end with eol character
#User must be sure this list is consistent with
#value of out_more
labels = ['sequence','nt','base_conf','puckering'] 

#open output file
fp = open(fname, "w")
#write dolumn labels
for l in labels:
    fp.write(l+delim)
fp.write(eol)
#begin program
for ind in range(len(json_data)):
    #find sequence
    f1=json_data[ind]['metadata']['str_id']+delim
    
    #determine number of nts and start loop
    for ind2 in dict_count(json_data[ind],'nts'):
        
        #f2 is index
        f2=find_it(json_data,ind,'nts',ind2,'index')+delim
        
        #f3 is baseSugar_conf                           
        f3=find_it(json_data,ind,'nts',ind2,'baseSugar_conf')+delim
        
        #f4 is puckering                                  
        f4=find_it(json_data,ind,'nts',ind2,'puckering')+delim
        f4s=f4
        f4l=f4


        if out_more == 1:
            fp.write(f1+f2+f3+f4l+eol)
            if out_screen == 1:
                print(f1+f2+f3+f4l)
        else:
            fp.write(f1+f2+f3+f4s+eol)
            if out_screen == 1:
                print(f1+f2+f3+f4s)
fp.close()

#code for json_detailedbp.py
 


#output file name
fname=".\\detailedbp.txt"
#following outputs all data to file when =1, 
#smaller output when =0
out_more=0 

#output to screen (=1) or silent (0)
out_screen=1

#field delimiter and eol
delim="#"
eol="\n"
 
#open json file  
json_fname = ".\\merged.json"    
with open(json_fname) as json_file:
    json_data = json.load(json_file) 

#names for column labels. should end with eol character
#User must be sure this list is consistent with
#value of out_more
labels = ['sequence','nt1','nt2','hbond1','hbond2','hbond3', 'hbond4'] 

#open output file
fp = open(fname, "w")
#write dolumn labels
for l in labels:
    fp.write(l+delim)
fp.write(eol)
    

results = []
for ind in range(len(json_data)):
        f1=json_data[ind]['metadata']['str_id']+delim
        
        for ind2 in dict_count(json_data[ind],'pairs'):
            f2=find_it(json_data,ind,'pairs',ind2,'nt1')
            f2s=f2[len(f2)-2:len(f2)-1]+delim 
            f2l=f2+delim+f2s
              
            f3=find_it(json_data,ind,'pairs',ind2,'nt2')
            f3s=f3[len(f3)-2:len(f3)-1]+delim
            f3l=f3+delim+f3s
                            
            f4=find_it(json_data,ind,'pairs',ind2,'hbonds_desc')+delim               
            #remove [  ] 
            f4=re.sub(r'\[.*?\]','', f4)
            #break at comma into two strings
            f4=re.sub(',', delim, f4)
            
            if out_more == 1:
                fp.write(f1+f2l+f3l+f4+eol)
                if out_screen == 1:
                    print(f1+f2l+f3l+f4)
            else:
                fp.write(f1+f2s+f3s+f4+eol)
                if out_screen == 1:
                    print(f1+f2s+f3s+f4)

fp.close()