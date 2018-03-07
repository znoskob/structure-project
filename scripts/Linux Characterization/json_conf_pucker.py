# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 04:37:55 2016

@author: kirkpacc
"""

def find_it (look_here,ind1,fld1,ind2,fld2):
    'finds value'
    if (fld1 in look_here[ind1]) and (fld2 in look_here[ind1][fld1][ind2]) :
        return str(look_here[ind1][fld1][ind2][fld2])
    else:
        return "blank#"   
        
def dict_count (dict,elem):
    'returns num of occurances of elem in dict'
    if elem in dict:
        return  range(len(dict[elem])) 
    else:
        return range(0,0);
        
#open json file 
import json    
json_fname = 'merged.json'        
with open(json_fname) as json_file:
    json_data = json.load(json_file)  
    
#following outputs all data to file when =1, 
#smaller output when =0
out_more=0 

#output to screen (=1) or silent (0)
out_screen=1

#output file name
fname="conf_pucker.txt"

#field delimiter and eol
delim="#"
eol="\n"

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
