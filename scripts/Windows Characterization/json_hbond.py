# -*- coding: utf-8 -*-
"""
Created on Mon Jul 25 04:37:55 2016

@author: kirkpacc
"""
    
import json
import re

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
 
#open json file  
json_fname = 'merged.json'        
with open(json_fname) as json_file:
    json_data = json.load(json_file) 

#following outputs all data to file when =1, 
#smaller output when =0
out_more=0 

#output to screen (=1) or silent (0)
out_screen=1

#output file name
fname="hbonds.txt"

#field delimiter and eol
delim="#"
eol="\n"

#names for column labels. should end with eol character
#User must be sure this list is consistent with
#value of out_more
labels = ['sequence','nt1','nt2','hbond1','hbond2','hbond3'] 

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
quit()

