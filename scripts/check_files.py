# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 09:17:09 2016
Modified on Thur Oct 20 11:37:05 2016 by kerichardson

@author: kirkpacc
"""

#finds files in directory that are missing bases
#assumes P, P01, P02 followed by nine primes (sugar)
#and then presence or absence of bases.  If total in 
#the section is >12, then assumed bases are present
#KR Also checks for multiple atom coordinates and
#KR files that only contain one line

def rough_check(g,content):
    #returns 0/1 for count is not/ok for group(g) in content
    count = 0
    for line in content:
        if (line[0:4] == "ATOM") and (int(line[24:27]) == g):
            count = count + 1
    if (count > 12) :
        result = 1
    else:
        result = 0        
    return result

def good_check(g,content,fname):
    #returns 0/1 for count is not/ok for group(g) in content
    count = 0
    review = []
    mult_atom = 0
    #KR Check for files with only one line
    if len(content) <= 1:
        is_empty = 1
        result = 0
        return result, is_empty
    else:
        is_empty = 0
    for line in content:
        #KR Check for multiple atom coordinates
        if ("A" or "B") in line[16:17]:
            mult_atom = mult_atom + 1
        if (line[0:4] == "ATOM" or line[0:6] == 'HETATM' or line[0:5] == 'HETAT') and (int(line[24:27]) == g):
            review.append(line)
            count = count + 1        
    #check first three lines for P
    #must check that values being tested are within range 
    p_ok = 0 
    tick_ok = 0        
    if count > 12:               
        for i in range(0,3):
            if "P" in review[i][:18]:
                p_ok = p_ok + 1
        #check next nine lines for '            
        for i in range(3,12):
            if "'" in review[i]:
                tick_ok = tick_ok +1
                #print(tick_ok)
            
    if ((count > 12) and (p_ok == 3) and (tick_ok == 9)and (mult_atom == 0) and (is_empty == 0)):
       result = 1
    else:
       result = 0
    return result,mult_atom,p_ok,tick_ok,count, is_empty

#KR Check for folder in path...if it does not exist, create it
def check_folder(folder):
    if os.path.isdir(folder) == False:
        print('Folder \'%s\' does not exist...creating' % folder)
        try:
            os.makedirs(folder)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise


import glob
import os
import shutil
import csv

script_path = ".\\"
clip_path = script_path+"clips\\" #relative path to clip directory
check_folder(script_path+'good_clips')
check_folder(script_path+'bad_clips')
good_files = script_path+"good_clips\\" #qssumes dir exists
bad_files = script_path+"bad_clips\\" #assumes dir exists
sequence = '*.pdb'
#sequence = '1MSY_CGUAAG_A.PDB'

#file_list = glob.glob(pdb_path)
file_list = glob.glob(clip_path+'%s'%sequence)

#loop through all files and test
split_line = []
bad_results = []
for fname in file_list:
    print("Checking "+fname)
    path,base = os.path.split(fname)
    file_ok = 1
    with open(fname) as f:
        content = f.readlines()
        if len(content) <= 1:
            is_empty = 1
            file_ok = 0
            bad_results.append([fname,0,0,0,0,0,0,1])
        else:
            is_empty = 0
            if (content[len(content)-3][24:27].strip()).isdigit():
                max_groups = int(content[len(content)-3][24:27])
                for g in range(1,max_groups):
                    group_ok,mult_atom,p_ok,tick_ok,count,is_empty = good_check(g,content,fname)
                    if group_ok == 0:
                        file_ok = 0
                        bad_results.append([fname,g,group_ok,p_ok,tick_ok,count,mult_atom,is_empty])
            else:
                max_groups = 0
                file_ok = 0
                bad_results.append([fname,0,0,0,0,0,0,0])
    if file_ok == 1: #write to good_files directory
        shutil.copy(fname,good_files+base)
        print("file ok. copy to good_clips "+fname)
    else:  #write to bad_files directgory  
        shutil.copy(fname,bad_files+base)  
        print("file not ok. copy to bad_clips "+fname)
#write badfile results to directory
with open(bad_files+"bad_file_data.txt", "w") as csv_file:
    writer = csv.writer(csv_file, delimiter='\t')
    writer.writerow(["File","group","file_ok","P Count","Sugar Count","Total Count","Total Mult Atom Count","IsEmpty"])
    for line in bad_results:
        writer.writerow(line)                           

                
           
            
        
    
