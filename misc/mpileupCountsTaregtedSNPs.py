#Brunilda Balliu
#March 16, 2016
# Make base count files from mpileup files

import sys

inFile = open(sys.argv[1],'r')

print 'CHR\tbp\tREF\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous\tDepth\tN_REF\t%_REF'
for line in inFile:
        data = line.strip().split('\t')
        if len(data)<5:
                line=inFile.next()
        elif data[3]==0:
               line=inFile.next() 
        else:

                CHR=data[0]
                bp = data[1]
                bases = data[4].upper()
                ref = data[2].upper()
        
                types = {'A':0,'G':0,'C':0,'T':0,'-':0,'+':[],'X':[]}

                i = 0
                while i < len(bases):
                        base = bases[i]
                        if base == '^' or base == '$':
                                i += 1
                        elif base == '-':
                                i += 1
                        elif base == '*':
                                types['-'] += 1
                        elif base == '+':
                                i += 1
                                addNum = int(bases[i])
                                addSeq = ''
                                for a in range(addNum):
                                        i += 1
                                        addSeq += bases[i]

                                types['+'].append(addSeq)
                        elif base == '.' or base == ',':
                                types[ref] += 1
                        else:
                                if types.has_key(base):
                                        types[base] += 1
                                else:
                                        types['X'].append(base)

                        i += 1

                adds = '.'
                if len(types['+']) > 0:
                        adds = ','.join(types['+'])

                amb = '.'
                if len(types['X']) > 0:
                        amb = ','.join(types['X'])

                depth=types['A']+types['G']+types['C']+types['T']+types['-']+len(types['+'])
                if depth<1:
                        line=inFile.next()
                        N_ref=types[ref]
                        prop_ref=0
                else:

                        N_ref=types[ref]
                        prop_ref=float(types[ref])/depth

                out = [CHR,bp,ref,types['A'],types['G'],types['C'],types['T'],types['-'],len(types['+']),adds,amb,depth,N_ref,prop_ref]
                print '\t'.join([str(x) for x in out])

