#!/usr/bin/env python
# coding: utf-8

import sys                      ##to use from command line, there are three variables
handle = open(sys.argv[1],'r')  #var 1 is the path to the HMMscan results file
out = sys.argv[2]               #var 2 is the path to the output file of unorganized data // var 3 is the output to the organized data
mydict, LRRlist, TIRlength, TIRlist = {}, [], 0, []
tircount, TIRpost = 0, 0
for line in handle:         ##records the first score necessary to determine the occurrence of a new sequence
    if line.startswith('#'):
        continue
    store = line.split()
    prevscore = float(store[7])
    break
for line in handle:
    if line.startswith('#'):
        continue
    store = line.split()   #stores line information in a list, every entry's information
    score = float(store[7])
    if score > prevscore:
        tircount, TIR_length, LRRlist, TIRpost, TIRlist = 0, 0, [], 0, [] ##reset containers and flags, new sequence
    prevscore = score
    if len(store) < 24 or len(store) > 25:  ##not LRR or TIR
        continue
    if len(store) == 24 and tircount == 1 and TIRpost == 1: #triggers if TIR domain previously written
        continue
    if len(store) == 24: ###triggers if TIR domain
        desc = store[22]+' '+store[23] ###desc gives identification of domain
        if desc != 'TIR domain':
            continue
        tircount = 1
        if store[9] != store[10]:
            key = store[3]
            TIR_length = int(store[18])-int(store[17])
            TIRlist.append(TIR_length)
            continue
        if store[9] == store[10]:
            key = store[3]
            TIR_length = int(store[18])-int(store[17])
            TIRlist.append(TIR_length)  ###calculate tir length, will only accept one domain, adds length of coding regions
            TIR_length = sum(TIRlist)/len(TIRlist)
            TIRpost = 1 #set TIRpost to 1, one TIR domain successfully recorded
            if key in mydict:
                mydict[key].append('TIR_Length:'+str(TIR_length)+'_'+str(store[10]))
                continue
            else:
                mydict[key] = ['TIR_Length:'+str(TIR_length)+'_'+str(store[10])]
                continue
    if len(store) == 25: ###triggers if LRR
        desc = store[22]+' '+store[23]+' '+store[24] ###desc gives identification of domain
        if desc != 'Leucine rich repeat':
            continue
        LRR_length = int(store[18])-int(store[17])
        LRRlist.append(LRR_length)
        key = store[3]
        if desc == 'Leucine rich repeat' and store[9] == store[10]:
            LRR_avgL = sum(LRRlist)/len(LRRlist)
            if key in mydict:
                mydict[key].append(str(LRR_avgL)+'_'+str(store[10]))
            else:
                mydict[key] = [str(LRR_avgL)+'_'+str(store[10])]
handle2 = open(out,'x')
for k,v in mydict.items():
    handle2.write(k+'\n')
    handle2.write(str(v)+'\n')
handle.close()
handle2.close()
################# Result Organization ########################
handle = open(out,'r')
handle2 = open(sys.argv[3]+'_Result_Verbose.txt','x')
TIR_lengths,LRR_lengths,LRR_counts = [],[],[]
count = 0
for line in handle:
    if 'TLR' in line:
        count += 1
        line = line.rstrip()
        handle2.write(line+'\n')
    else:
        count += 1
        line = line.replace('[','')
        line = line.replace(']','')
        line = line.replace(',','')
        line = line.replace("'",'')
        store = line.split()
        for ele in store:
            if 'TIR' in ele:
                num_s = ele.index(':')+1
                TIR_l_n = ele[num_s:]
                TIR_lengths.append(float(TIR_l_n[:TIR_l_n.index('_')]))
            elif 'TIR' not in ele:
                div = ele.index('_')
                LRR_l = ele[:div]
                LRR_lengths.append(float(LRR_l))
                LRR_r = ele[div+1:]
                LRR_counts.append(int(LRR_r))
            else:
                continue
        TIR_l_avg = sum(TIR_lengths)/len(TIR_lengths)
        if LRR_lengths != []:
            LRR_l_avg = sum(LRR_lengths)/len(LRR_lengths)
        if LRR_counts != []:
            LRR_c_avg = sum(LRR_counts)/len(LRR_counts)
        if LRR_counts == [] and LRR_lengths == []:
            handle2.write('TIR Length AVG: '+str(TIR_l_avg)+'\n')
            TIR_lengths = []
            continue
        handle2.write('TIR Length AVG: '+str(TIR_l_avg)+' '+'LRR Count AVG: '+str(LRR_c_avg)+' '+'LRR Length AVG: '+str(LRR_l_avg)+'\n')
        TIR_lengths,LRR_lengths,LRR_counts = [],[],[]
handle.close()
handle2.close()
handle3 = open(sys.argv[3]+'_Result_Verbose.txt','r')   ##organizes data as such: name, TIR_Length, LRR_Count, LRR_Avg_Length
handle4 = open(sys.argv[3]+'.csv','x')
for line in handle3:
    line = line.rstrip()
    if 'TLR' in line:
        query = line
        continue
    store = line.split()
    if len(store) == 4:
        bundle = store[3]
        handle4.write(query+','+bundle+'\n')
    if len(store) == 12:
        bundle = store[3]+','+store[7]+','+store[11]
        handle4.write(query+','+bundle+'\n')
handle3.close()
handle4.close()