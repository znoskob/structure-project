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
        return "blank"   
        
def dict_count (dict,elem):
    'returns num of occurances of elem in dict'
    if elem in dict:
        return  range(len(dict[elem])) 
    else:
        return range(0,0);
        
def find_in_str(str):
    #returns the single occurence of a string in list that is in str
    #exits if none of more than one string in list is found in str 
    #list is given below
    to_look_for=['<<','<>','><','>>']
    x=[int(str.find(to_look_for[0])),int(str.find(to_look_for[1])),int(str.find(to_look_for[2])),int(str.find(to_look_for[3]))]
    #print (x)
    if x == [-1,-1,-1,-1]: 
        #not found
        return "";
    else:
        return to_look_for[x.index(max(x))]
        
#    if max((x))-3 == sum((x)) :
#        pos=x.index(max(x)); 
#        return to_look_for[pos]
#    else:
#        print ('not ok in find_in_str')
#        exit()

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

#output not included
out_not_included = 1

#output file name
fname="stacks.txt"

#field delimiter and eol
delim="#"
eol="\n"

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
    #print("ind: "+str(ind))
    #print("json_data[ind]: "+str(json_data[ind]))
    f1=json_data[ind]['metadata']['str_id']+delim
    
    #determine number of nonPairs. Loop to find nt1,nt2,stacking
    for ind2 in dict_count(json_data[ind],'nonPairs'):
 
        #f2 is nt1
        f2=find_it(json_data,ind,'nonPairs',ind2,'nt1')
        f2s = f2[len(f2)-2:len(f2)-1]+delim   
        f2l=f2 + f2[len(f2)-2:len(f2)-1]+delim
        
        #f3 is nt2
        f3=find_it(json_data,ind,'nonPairs',ind2,'nt2')
        f3s = f3[len(f3)-2:len(f3)-1]+delim 
        f3l=f3 + f3[len(f3)-2:len(f3)-1]+delim
        
        #f4 is stacking
        f4=find_it(json_data,ind,'nonPairs',ind2,'stacking')+delim
        f4s = find_in_str(f4) 
        f4l = f4 + f4s +delim
        
        outstringl=f1+f2+f3+f4
        
        #output to file only if all fields are not blank
        #output to file only if F4s != ""
        if f4s != "":
            if out_more == 1:
                fp.write(f1+f2l+f3l+f4l+eol)
                if out_screen == 1:
                    print(f1+f2l+f3l+f4l)
            else:
                fp.write(f1+f2s+f3s+f4s+eol)
                if out_screen == 1:
                    print(f1+f2s+f3s+f4s) 
        elif out_not_included == 1: 
            if out_more == 1:
                fp.write(f1+f2l+f3l+f4l+"not included"+eol)
                if out_screen == 1:
                    print(f1+f2l+f3l+f4l+"not included") 
            else:
                fp.write(f1+f2s+f3s+f4s+"not included"+eol)
                if out_screen == 1:  
                    print(f1+f2s+f3s+f4s+"not included")                    
fp.close()

        


