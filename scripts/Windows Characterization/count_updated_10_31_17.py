from __future__ import division
from operator import itemgetter
import csv
import pandas as pd
from collections import Counter
import itertools
from tabulate import tabulate
import os
import re
import datetime



cwd = os.path.split(os.getcwd())[1]
cutoff = os.path.split(os.getcwd())[0]
cutoff1 = os.path.basename(os.path.normpath(cutoff))
data_list = re.findall(r'\d+(?:\.\d+)?', cutoff1)
cutoff2 = str(data_list).strip('[]')
cutoff3 = cutoff2.replace("'", "")

now = datetime.datetime.now()


d = pd.read_csv('stacks.csv')
unique_sequences = list(d['sequence'].unique())

total = len(unique_sequences)

print '\n'
print "Date and time of run:", now.strftime("%Y-%m-%d %H:%M")
print 'Degenerate Seqence:', cwd
print 'Total Number of Structures:', total
print  'RMSD Cutoff:', cutoff3,'A'
print '\n'

nt1a = list(d['nt1'])
ntlength = len(nt1a)

if ntlength == 0:
    print 'No stacking interactions found.'

if ntlength != 0:
    d['nts_stacking'] = d[['nt1','nt2','stacking']].apply(lambda x : '{}{}{}'.format(x[0],x[1],x[2]), axis=1)
    stacks = list(d['nts_stacking'])
    P = Counter(stacks)
    
    stacksab = list(d['nts_stacking'])
    stacksc = [x for x in stacksab if not 'not included' in x]
    
    P = Counter(stacksc)
    for x in P:
        P[x]/=total

    for x in P:
        P[x]*=100

    headers = ['NT 1, NT 2 & Stacking', 'Percentages']
    data = sorted([(k,v) for k,v in P.iteritems() if v != 0])
    print(tabulate(data, headers=headers, tablefmt='grid'))

    da = pd.DataFrame.from_dict(P, orient='index', dtype=Counter)
    da.to_csv('stacks_percent.csv')
    dab = pd.read_csv('stacks_percent.csv', names = ["NT 1, NT 2 & Stacking", "Percentage of Structures with Stacking Interaction"])

print "\n"




df = pd.read_csv('conf_pucker.csv')

nt = list(df['nt'])

df['nt_conf'] = df[['nt', 'base_conf']].apply(lambda x : '{}{}'.format(x[0],x[1]), axis=1)

conf = list(df['nt_conf'])

conf1 = [x for x in conf if not 'nan' in x]
conf2 = [x for x in conf1 if not 'blank' in x]



ntmax = max(nt)

while ntmax > 0:
    
    nt_count = nt.count(ntmax)
    conf_count = Counter(conf2)

    ntmax = ntmax - 1



    for x in conf_count:
        conf_count[x] = conf_count[x]/nt_count
        conf_count[x] = conf_count[x]*100



headers1 = ['NT & Base Conformation', 'Percentages']
data = sorted([(k,v) for k,v in conf_count.iteritems() if v != 0])
print(tabulate(data, headers=headers1, tablefmt='grid'))
print '\n'




dir_path = os.path.dirname(os.path.realpath(__file__))

df['nt_pucker'] = df[['nt','puckering']].apply(lambda x : '{}{}'.format(x[0],x[1]), axis=1)

pucker = list(df['nt_pucker'])

pucker1 = [x for x in pucker if not 'nan' in x]
pucker2 = [x for x in pucker1 if not 'blank' in x]

ntmax1 = max(nt)

while ntmax1 > 0:
    nt_count1 = nt.count(ntmax1)
    pucker_count = Counter(pucker2)

    ntmax1 = ntmax1 - 1

    for x in pucker_count:
        pucker_count[x] = pucker_count[x]/nt_count1
        pucker_count[x] = pucker_count[x]*100

headers2 = ['NT & Sugar Pucker', 'Percentages']
data = sorted([(k,v) for k,v in pucker_count.iteritems() if v != 0])
print(tabulate(data, headers=headers2, tablefmt='grid'))
print '\n'


dfb = pd.DataFrame.from_dict(conf_count, orient='index', dtype=Counter)
dfb.to_csv('conf_percent.csv')
dfbb = pd.read_csv('conf_percent.csv', names = ["NT & Base Conformation", "Percentage of NTs with Conformation"])

dfa = pd.DataFrame.from_dict(pucker_count, orient='index', dtype=Counter)
dfa.to_csv('pucker_percent.csv')
dfab = pd.read_csv('pucker_percent.csv', names = ["NT & Sugar Pucker", "Percentage of NTs with Pucker"])


df1 = pd.read_csv('hbonds.csv',error_bad_lines=False)

sequences = list(df1['sequence'].unique())

total = len(sequences)

nts_hbond1a = list(df1['hbond1'])
totalhbonds = len(nts_hbond1a)

if totalhbonds == 0:
    print 'No non-pairing hydrogen bonding interactions found.'

if totalhbonds != 0:
    hbond1 = {k: v for v, k in enumerate(nts_hbond1a)}
    df1['nts_hbond1'] = df1[['nt1','nt2','hbond1']].apply(lambda x : '{}{}{}'.format(x[0],x[1],x[2]), axis=1)
    df1['nts_hbond2'] = df1[['nt1','nt2','hbond2']].apply(lambda x : '{}{}{}'.format(x[0],x[1],x[2]), axis=1)
    df1['nts_hbond3'] = df1[['nt1','nt2','hbond3']].apply(lambda x : '{}{}{}'.format(x[0],x[1],x[2]), axis=1)


    hbonds0 = list(df1['nts_hbond1'])
    hbonds1 = list(df1['nts_hbond2'])
    hbonds2 = list(df1['nts_hbond3'])

    hbond = itertools.chain(hbonds0, hbonds1, hbonds2)
    hbond_list = list(hbond)
    

    hbonds = [x for x in hbond_list if not 'nan' in x]
    hbonds1a = [x for x in hbonds if not 'blank' in x]
    count16 = Counter(hbonds1a)

    for x in count16:
        count16[x]/=total
    
    for x in count16:
        count16[x]*=100

    headers3 = ['NT 1, NT 2, and Hydrogen Bond', 'Percentages']
    data = sorted([(k,v) for k,v in count16.iteritems() if v !=0])
    print(tabulate(data, headers=headers3, tablefmt='grid'))


    df1a = pd.DataFrame.from_dict(count16, orient='index', dtype=Counter)
    df1a.to_csv('hbond_percent.csv')
    df1ab = pd.read_csv('hbond_percent.csv', names = ["NT 1, NT 2 & Hydrogen Bond", "Percentage of Structures with Non-Pairing HB"])
print '\n'


df3 = pd.read_csv('bp.csv')
sequences_bp = list(df3['sequence'].unique())

totals = len(sequences_bp)

nt1aa = list(df3['nt1'])

totalbp = len(nt1aa)

if totalbp == 0:
    print 'No base pairing interactions found.'
if totalbp != 0:
    df3['nts_bp'] = df3[['nt1','nt2']].apply(lambda x : '{}{}'.format(x[0],x[1]), axis=1)

    pairs = list(df3['nts_bp'])

    pairs1 = [x for x in pairs if not 'nan' in x]
    pairs2 = [x for x in pairs1 if not 'blank' in x]

    counta = Counter(pairs2)


    for x in counta:
        counta[x]/=totals
        counta[x]*=100


    headers4 = ['Paired Nucleotides', 'Percentages']
    data = sorted([(k,v) for k,v in counta.iteritems() if v != 0])
    print(tabulate(data, headers=headers4, tablefmt='grid'))
    print '\n'


    df3['nts_bp_type'] = df3[['nt1','nt2','name']].apply(lambda x : '{}{}{}'.format(x[0],x[1],x[2]), axis=1)

    pair = list(df3['nts_bp_type'])

    pair1 = [x for x in pair if not 'nan' in x]
    pair2 = [x for x in pair1 if not 'blank' in x]

    count1 = Counter(pair2)


    for x in count1:
        count1[x]/=totals
        count1[x]*=100


    headers5 = ['NT 1, NT 2 & Base Pair Classification', 'Percentages']
    data = sorted([(k,v) for k,v in count1.iteritems() if v != 0])
    print(tabulate(data, headers=headers5, tablefmt='grid'))

    df3b = pd.DataFrame.from_dict(counta, orient='index', dtype=Counter)
    df3b.to_csv('bp_percent.csv')
    df3bb = pd.read_csv('bp_percent.csv', names = ["Paired Nucleotides", "Percentage of Structures with Pairing"])


    df3a = pd.DataFrame.from_dict(count1, orient='index', dtype=Counter)
    df3a.to_csv('bpclass_percent.csv')
    df3ab = pd.read_csv('bpclass_percent.csv', names = ["NT 1, NT 2 & Base Pair Classification", "Percentage of Base Pair Classification"])
print '\n'




df4 = pd.read_csv('detailedbp.csv',error_bad_lines=False)

sequences = list(df1['sequence'].unique())

totalb = len(sequences)
hbond1aa = list(df4['nt1'])
totalbphb = len(hbond1aa)

if totalbphb != 0:
    hbond1 = {k: v for v, k in enumerate(hbond1aa)}
    df4['nts_hbond1'] = df4[['nt1','nt2','hbond1']].apply(lambda x : '{}{}{}'.format(x[0],x[1],x[2]), axis=1)
    df4['nts_hbond2'] = df4[['nt1','nt2','hbond2']].apply(lambda x : '{}{}{}'.format(x[0],x[1],x[2]), axis=1)
    df4['nts_hbond3'] = df4[['nt1','nt2','hbond3']].apply(lambda x : '{}{}{}'.format(x[0],x[1],x[2]), axis=1)
    df4['nts_hbond4'] = df4[['nt1','nt2','hbond4']].apply(lambda x : '{}{}{}'.format(x[0],x[1],x[2]), axis=1)
    keywordFilter = set(['blank'])

    hbonds3 = list(df4['nts_hbond1'])
    hbonds4 = list(df4['nts_hbond2'])
    hbonds5 = list(df4['nts_hbond3'])
    hbonds6 = list(df4['nts_hbond4'])
    

    hbond0 = itertools.chain(hbonds3, hbonds4, hbonds5, hbonds6)
    hbond_lists = list(hbond0)



    hbonds1 = [x for x in hbond_lists if not 'nan' in x]
    hbonds2 = [x for x in hbonds1 if not 'blank' in x]

    countb = Counter(hbonds2)
    
    for x in countb:
        countb[x]/=totalb

    for x in countb:
        countb[x]*=100

    headers6 = ['NT 1, NT 2 & Pairing Hydrogen Bonds', 'Percentages']
    data = sorted([(k,v) for k,v in countb.iteritems() if v != 0])
    print(tabulate(data, headers=headers6, tablefmt='grid'))

    df4a = pd.DataFrame.from_dict(countb, orient='index', dtype=Counter)
    df4a.to_csv('phbond_percent.csv')
    df4ab = pd.read_csv('phbond_percent.csv', names = ["NT 1, NT 2 & Pairing Hydrogen Bonds", "Percentage of Pairing Hydrogen Bonds"])
print '\n'

if totalbp != 0:
    total_frames = pd.concat([dab, dfbb, dfab, df1ab, df3bb, df3ab, df4ab], axis=1)
    total_frames.to_csv('summary1.csv', index=None)

    first_row_num = 1
    delete_row = {2}
    with open('summary1.csv', 'rt') as infile, open('summary.csv', 'wt') as outfile:
        outfile.writelines(row for row_num, row in enumerate(infile, first_row_num)
                       if row_num not in delete_row)




