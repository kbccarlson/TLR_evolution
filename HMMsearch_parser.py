import sys									### to use in command line, the first argument is the path to the HMMsearch output file
                                            ### the second argument is the path of the source of sequences searched with HMMsearch, e.x. 'Coquerels_sifaka.faa'
from Bio import SearchIO as sio 			###use searchIO to format data

x = sio.read(sys.argv[1],format='hmmer3-text')
f = open(sys.argv[1]+'_summary.txt','x')					###create a summary file to dump info
flag = 0
for line in x.hsps:				###.hsps is an attribute that contains the information we want in the hit table
	line = str(line)
	origin = line.find('evalue')	###find e value for origin point reference
	if 'e' not in line[origin+6:origin+16]:  
		continue
	e = line[origin+6:].find('e')   ###Find e, relevant to origin
	e = (origin+6) + e ###add origins pointer to the previous e location derived from origin, e now points to e without variable reference
	prefix = float((line[origin+7:e]))  
	suffix = int((line[e+1:e+4]))		
	evalue = (prefix)*(10**(suffix))	###calc e value in meaningful way
	if evalue > .01: 					###inclusion threshhold filtering
		continue
	if flag == 0:						
		f.write(line)
		f.write('\n')
		flag = 1
		continue
	if flag != 0:
		m = line.find('Query:')			
		line = line[:m]+'\n'+'     '+line[m:]   
		f.write(line)
		f.write('\n')
		continue
f.close()

handlef = open(sys.argv[1]+'_summary.txt','r')
mydict, myhits = {}, []
for line in handlef:
    line = line.rstrip()
    if 'Hit:' in line:
        start = line.index(':')
        if '.1' in line:
            end = line.index('.1')
        if '.2' in line:
            end = line.index('.2')
        hit = line[(start+2):(end+3)]
        myhits.append(hit)
mydict['record']=myhits
y = list(mydict.items())
handlef.close()

fasta_list = [] ##holds sequence and header info
flag = 0
for q in range(len(y[0][1])):
    f = open(sys.argv[2],"r")
    pattern = str(y[0][1][q]) 
    for line in f:              ###append matches with their sequences
        if pattern in line:
            fasta_list.append(line)
            flag = 1
            continue
        if flag == 1 and '*' in line:  
            fasta_list.append(line)
            flag = 0
            continue            
        if flag == 1 and '>' in line and pattern not in line: 
            flag = 0
            continue
        if flag == 1:
            fasta_list.append(line)
            continue
        else:
            continue
    f.close()
                
k = open(sys.argv[1]+"_hits.txt","x")   
flag2 = 0      
for j in fasta_list:
    j = j.replace("\n","")  
    if '>' in j and flag2 == 0:  
        k.write(j+'\n')
        flag2 = 1
    elif '>' in j and flag2 != 0: 
        k.write('\n'+j+'\n')
    elif '>' not in j: 
        k.write(j)
k.close()