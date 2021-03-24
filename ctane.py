"""------------------------------------------------------------------------------------------
TANE Algorithm for discovery of exact conditional functional dependencies
Author: Nabiha Asghar, nasghar@uwaterloo.ca
March 2015
Use for research purposes only.
Please do not re-distribute without written permission from the author
Any commerical uses strictly forbidden.
Code is provided without any guarantees.
----------------------------------------------------------------------------------------------"""

def replace_element_in_tuple(tup, elementindex, elementval):
    if type(elementval)==tuple:
        elementval = elementval[0]
    newtup = list(tup)
    newtup[elementindex] = elementval
    newtup = tuple(newtup)
    return newtup

def add_element_in_tuple(spxminusa, ca):
    thelist = list(spxminusa)
    thelist.append(ca[0])
    return tuple(thelist)

def validcfd(xminusa, x, a, spxminusa, sp, ca):
    global dictpartitions
    if xminusa == '' or a == '': 
        return False
    indexofa = x.index(a)
    newsp0 = add_element_in_tuple(spxminusa, ca)
    newsp1 = replace_element_in_tuple(sp, indexofa, ca)   #this is sp, except that in place of value of a we put ca
    if (x, newsp1) in dictpartitions.keys():
        if  len(dictpartitions[(xminusa, spxminusa)]) == len(dictpartitions[(x, newsp1)]):# and twodlen(dictpartitions[(xminusa, spxminusa)]) == twodlen(dictpartitions[(x, newsp1)]):
            return True    
    return False

def twodlen(listoflists):
	summ = 0
	for item in listoflists:
		summ = summ + len(item)
	return summ

def greaterthanorequalto(upxminusa, spxminusa): # this is actually greaterthan or equal to
    if upxminusa == spxminusa: 
        return True
    flag = True
    for index in range(0, len(upxminusa)):
        if not (spxminusa[index]=='--'):
            if (not (upxminusa[index] == spxminusa[index])):
                flag = False
    return flag

def doublegreaterthan(upxminusa, spxminusa): 
    if upxminusa == spxminusa: 
        return False
    flag = True
    for index in range(0, len(upxminusa)):
        if (not spxminusa[index]=='--'):
            if (not (upxminusa[index] == spxminusa[index])):
                flag = False
    return flag
    
def compute_dependencies(level, listofcols):
    global dictCplus
    global finallistofCFDs
    global listofcolumns
    for (x,sp) in level:
        for a in x:
            for (att, ca) in dictCplus[(x, sp)]:
                if att == a:
                    newtup =  spXminusA(sp, x, a)      ### tuple(y for y in sp if not sp.index(y)==x.index(a)) # this is sp[X\A]                             
                    if validcfd( x.replace(a,''), x, a, newtup, sp, ca) and not ([x.replace(a,''), a, [newtup, ca]] in finallistofCFDs):
                        finallistofCFDs.append([x.replace(a,''), a, [newtup, ca]])
                        for (xx, up) in level:
                            if xx==x:
                                newtup0 =  spXminusA(up, x, a)          ### tuple(y for y in up if not up.index(y)==x.index(a)) # this is up[X\A]
                                if up[x.index(a)]==ca[0] and greaterthanorequalto(newtup0, newtup) :
                                    if (a, ca) in dictCplus[(x,up)]: dictCplus[(x,up)].remove((a,ca))
                                    listofcolscopy = listofcols[:]
                                    for j in x: # this loop computes R\X
                                        if j in listofcolscopy: listofcolscopy.remove(j)
                                    for b_att in listofcolscopy: # this loop removes each b in R\X from C+(X,up)
                                        stufftobedeleted = []
                                        for (bbval, sometup) in dictCplus[(x,up)]:
                                            if b_att == bbval:
                                                stufftobedeleted.append((bbval,sometup))                        
                                        for item in stufftobedeleted:
                                            dictCplus[(x,up)].remove(item)

def prune(level):
    global dictCplus
    stufftobedeleted=[]
    for (x,sp) in level:
        if len(dictCplus[(x,sp)])==0:
            stufftobedeleted.append((x,sp))
    for item in stufftobedeleted:
        level.remove(item)

def computeCplus(level): # for each tuple (x,sp) in the list level, it computes C+(x,sp), which is a list of (attribute, value) tuples) 
    global listofcolumns
    global dictCplus
    listofcols = listofcolumns[:]
    for (x,sp) in level: #sp is a tuple of strings like this: ('aa', 'bb', 'cc') or ('aa', )     
        thesets=[]
        for b in x:
            indx = x.index(b) # the index where b is located in x
            spcopy =  spXminusA(sp, x, b)     ### tuple(y for y in sp if not sp.index(y)==indx)
            spcopy2 = sp[:]            
            if (x.replace(b,''), spcopy ) in dictCplus.keys():
                temp = dictCplus[(x.replace(b,''), spcopy)]
            else: temp = []   # is this correct???? should I put [] here?
            thesets.insert(0, set(temp))
        if list(set.intersection(*thesets)) == []:
            dictCplus[(x,sp)] = []
        else:
            dictCplus[(x,sp)] = list(set.intersection(*thesets))

def initial_Cplus(level):
    global listofcolumns
    global dictCplus
    computeCplus(level)
    for (a,ca) in level:
        stufftobedeleted = []
        for (att, val) in dictCplus[(a,ca)]:
            if att==a and not val==ca:
                stufftobedeleted.append((att,val))
        for item in stufftobedeleted:
            dictCplus[(a,ca)].remove(item)

def populateL1(listofcols):    
    global k_suppthreshold
    l1 = []
    attributepartitions = computeAttributePartitions(listofcols)
    for a in listofcols:
        l1.append((a, ('--',)))
        for eqclass in attributepartitions[a]:
            if len(eqclass)>= k_suppthreshold:
                l1.append( (a, (str(data2D.iloc[eqclass[0]][a]) , ) ) )
    computeInitialPartitions(l1, attributepartitions) # populates the dictpartitions with the initial partitions (X,sp) where X is a single attribute
    return l1

def computeInitialPartitions(level1, attributepartitions):
	global data2D
	global dictpartitions # dictpartitions[(x,sp)] is of the form [[0,1,2]]. So simply a list of lists of indices  
	for (a,sp) in level1:
		dictpartitions[(a,sp)]=[]
		dictpartitions[(a,sp)] = attributepartitions[a]

def old_computeInitialPartitions(level1, attributepartitions):
    global data2D
    global dictpartitions # dictpartitions[(x,sp)] is of the form [[0,1,2]]. So simply a list of lists of indices  
    for (a,sp) in level1:
        dictpartitions[(a,sp)]=[]
        if sp[0]=='--':
            dictpartitions[(a,sp)] = attributepartitions[a]
        else:
            for eqclass in attributepartitions[a]:
                if str(data2D.iloc[eqclass[0]][a])==sp[0]:
                    dictpartitions[(a,sp)].append(eqclass)

def computeAttributePartitions(listofcols): # compute partitions for every attribute 
    global data2D    
    attributepartitions = {}
    for a in listofcols:
        attributepartitions[a]=[]
        for element in list_duplicates(list(data2D[a])): # list_duplicates returns 2-tuples, where 1st is a value, and 2nd is a list of indices where that value occurs
            if len(element[1])>0: # if >1, then ignore singleton equivalence classes
                attributepartitions[a].append(element[1])
    return attributepartitions

def list_duplicates(seq):
    tally = defaultdict(list)
    for i,item in enumerate(seq):
        tally[item].append(i)
    return ((key,locs) for key,locs in tally.items() 
                            if len(locs)>0)

def sometuplematchesZUP(z,up):
    global dictpartitions
    global k_suppthreshold
    sumofmatches = 0
    for eqclass in dictpartitions[(z, up)]:
        sumofmatches = sumofmatches +  len(eqclass)
    if sumofmatches >= k_suppthreshold:
        return True
    else:
        return False

def generate_next_level(level):
    nextlevel=[]
    for i in range(0,len(level)): # pick an element
        for j in range(i+1, len(level)): # compare it to every element that comes after it. 
            if ((not level[i][0]==level[j][0]) and level[i][0][0:-1]==level[j][0][0:-1] and level[i][1][0:-1]==level[j][1][0:-1]):
                z = level[i][0] + level[j][0][-1]
                up = tuple(list(level[i][1]) + [level[j][1][-1]])
                (z, up) = sortspbasedonx(z, up)
                partition_product((z,up), level[i], level[j])
                if sometuplematchesZUP(z,up):
                    flag = True
                    for att in z:
                        indexofatt = z.index(att) # where is att located in z                        
                        up_zminusa = spXminusA(up, z, att)
                        zminusa = z.replace(att,'')
                        if not ((zminusa, up_zminusa) in level):
                            flag = False
                    if flag:
                        nextlevel.append((z, up))
    return nextlevel

def spXminusA(sp, x, a):
    indexofa = x.index(a)
    mylist=[]
    for i in range(0, len(sp)):
        if not i==indexofa:
            mylist.append(sp[i])
    return tuple(mylist)

def partition_product(zup, xsp, ytp):
    global dictpartitions
    global tableT
    tableS = ['']*len(tableT)
    partitionXSP = dictpartitions[xsp]
    partitionYTP = dictpartitions[ytp]
    partitionZUP = []
    for i in range(len(partitionXSP)):
        for t in partitionXSP[i]:
            tableT[t] = i
        tableS[i]=''
    for i in range(len(partitionYTP)):
        for t in partitionYTP[i]: 
            if ( not (tableT[t] == 'NULL')): 
                tableS[tableT[t]] = sorted(list(set(tableS[tableT[t]]) | set([t]))) 
        for t in partitionYTP[i]: 
            if (not (tableT[t] == 'NULL')) and len(tableS[tableT[t]])>= 1 : 
                partitionZUP.append(tableS[tableT[t]]) 
            if not (tableT[t] == 'NULL'): tableS[tableT[t]]='' 
    for i in range(len(partitionXSP)): 
        for t in partitionXSP[i]: 
            tableT[t]='NULL'
    dictpartitions[zup] = partitionZUP
    dictpartitions[zup] = partitionZUP

def sortspbasedonx(x,sp):
    x = list(x)
    points = zip(x,sp)
    sorted_points = sorted(points)
    new_x = [point[0] for point in sorted_points]
    new_sp = [point[1] for point in sorted_points]
    return (''.join(new_x), tuple(new_sp))

def correlation_matrix(data_path):   
    data = pandas.read_csv(data_path) 
    x = data.phik_matrix()
    x = x.replace(to_replace = 1, value = 0) # this removes any 1:1 corrolation ie: age:age, experience:experience etc.
    return x

def remove_redundant_corrs(corr_matrix):
    """ Returns a dict of all the feature pair correlations """
    corr_dict = corr_matrix.unstack().to_dict()

    filtered_corrs = {}
    for col_pair, corr in list(corr_dict.items()):
        if corr == 0 or \
        sorted(col_pair) in [sorted(key) for key in list(filtered_corrs.keys())]: 
            continue

        filtered_corrs[col_pair] = corr

    return filtered_corrs
#------------------------------------------------------- START ---------------------------------------------------

if __name__ == "__main__":
    from collections import defaultdict
    import numpy as NP
    import pandas
    import sys
    import os
    import phik

    os.environ["MODIN_ENGINE"] = "dask"

    from modin.pandas import *
    # for me it's win32. if you are using linux, please try using:
    # import multiprocessing.popen_spawn_posix
    import multiprocessing.popen_spawn_win32
    from distributed import Client
    client = Client()


    if len(sys.argv) > 1:
        infile=str(sys.argv[1])
    # if len(sys.argv) > 2:
    #     k=int(sys.argv[2])

    data2D = read_csv(infile)

    thresh = lambda corr: 2 if corr >=0.46 else 2

    corr_matrix = correlation_matrix(infile)
    filtered_corrs = remove_redundant_corrs(corr_matrix)
    mapped_col_k_vals = {col_pair: thresh(corr) for col_pair, corr in filtered_corrs.items()}
    
    for col_pair, k in mapped_col_k_vals.items():
        totaltuples = len(data2D.index)
        listofcolumns = list(col_pair) # returns ['A', 'B', 'C', 'D', .....]
        tableT = ['NULL']*totaltuples # this is for the table T used in the function partition_product
        k_suppthreshold = k
        L0 = []

        dictpartitions = {} # maps 'stringslikethis' to a list of lists, each of which contains indices
        finallistofCFDs=[]
        L1=populateL1(listofcolumns[:])  # L1 is a list of tuples of the form [ ('A', ('val1') ), ('A', ('val2') ), ..., ('B', ('val3') ), ......]
        dictCplus = {('',()): L1[:]}
        l=1
        L = [L0,L1]

        while (not (L[l] == [])):
            if l==1:
                # print("l == 1")
                initial_Cplus(L[l])
            else:
                # print("l != 1")
                computeCplus(L[l])
            # print("computing the dependencies")
            compute_dependencies(L[l],listofcolumns[:])
            # print("pruning..")
            prune(L[l])
            # print("generating the next level")
            temp = generate_next_level(L[l])
            L.append(temp)
            l=l+1
            # print("One iteration done...\n")       

        # have to manually map columns...
        column_mapping = {
            "A": "Age",
            "B": "Workclass",
            "C": "fnlwgt",
            "D": "Education",
            "E": "Education-num",
            "F": "Marital Status",
            "G": "Occupation",
            "H": "Relationship",
            "I": "Race",
            "J": "Gender",
            "K": "Capital-Gain",
            "L": "Capital-Loss",
            "M": "HoursPerWeek",
            "N": "Native-Country",
            "O": "Class"
        }

        for CFD in finallistofCFDs:
            for i, cols in enumerate(CFD[:2]):
                mapping = ""
                for col in cols:
                    mapping += f"{column_mapping[col]}, "
                CFD[i] = mapping.rstrip(", ")

        print(f"List of {col_pair}'s CFDs: " , finallistofCFDs)
        print("Total number of CFDs found: ", len(finallistofCFDs))
