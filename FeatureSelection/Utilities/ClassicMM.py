from multiprocessing import Pool, Semaphore
import time


def init(s):
    global semaphore
    semaphore = s
    return

def ExistDseparator(TargetNode,Xi, Z, X, alpha):
    flagExist = False
    dsepSet=[]
    #counter=0
    #print_names(Z)
    for i in range(1,(2**len(Z))):
        IDsubsetZ_dec = i
        IDsubsetZ_bin = bin(IDsubsetZ_dec)
        subsetZ = getZsubset(IDsubsetZ_bin,Z)
        # no cache
        #print("from exist dseparator")
        dep = Dep(TargetNode,Xi,subsetZ, X, alpha)
        #print(subsetZ)
        #print(dep)
        if (dep==0):
            flagExist = True
            dsepSet = subsetZ
            break
    #print("Module exist d-separator: ",counter)
    return [flagExist,dsepSet]


# def Dep_batch(params_arr):
#     returnable = []
#     for e in params_arr:
#         if (s.get_value()==1):
#             result = Dep(*e)
#             if result == 0:
#                 l.acquire()
#                 s.acquire()
#                 time.sleep(0.5)
#                 s.release()
#                 l.release()
#                 return [True,e[2]]
#         else:
#             return [False,[]]

def ExistDseparator_concurrent(TargetNode,Xi, Z, X, alpha):
    flagExist = False
    dsepSet=[]
    #counter=0
    #print_names(Z)
    if (len(Z)>10):
        n_processes = 8
        s = Semaphore(1)
        p = Pool(n_processes,initializer=init,initargs=(s))
        params_arr = [[TargetNode,Xi,getZsubset(bin(index),Z), X, alpha] for x in list(range(1,2**len(Z)))]
        batch_arr = []
        step = 2**len(Z)/n_processes
        for n in range(0,len(n_processes)):
            batch_arr.append(params_arr[n:n+step])
        results = p.map(Dep_batch,params_arr)
        for r in results:
            if r[0] == True:
                return r
        return [False, []]
        
    for i in range(1,(2**len(Z))):
        IDsubsetZ_dec = i
        IDsubsetZ_bin = bin(IDsubsetZ_dec)
        subsetZ = getZsubset(IDsubsetZ_bin,Z)
        # no cache
        #print("from exist dseparator")
        dep = Dep(TargetNode,Xi,subsetZ, X, alpha)
        #print(subsetZ)
        #print(dep)
        if (dep==0):
            flagExist = True
            dsepSet = subsetZ
            break
    #print("Module exist d-separator: ",counter)
    return [flagExist,dsepSet]

def getZsubset(bin_id,Z):
    bin_str=str(bin_id[::-1])
    Zsubset=[]
    for i in range(0,len(bin_str)):
        if bin_str[i]=='1':
            Zsubset.append(Z[i])
    return Zsubset

def MaxMinHeuristic(TargetNode, CPC, Universe, X, alpha):
    F=[]
    assocF=-1
    Z = CPC
    fixedCondVars = []
    if (len(CPC)>0):
        Z = CPC[0:-1]           # all but the last one   
        fixedCondVars = [CPC[-1]] # we use last one
    for i in range(len(Universe)-1,-1,-1):
        if (len(Universe[i])==0):
            continue
        assoc = MinAssoc(TargetNode,Universe[i],Z,fixedCondVars,X,alpha)
        if (assoc>assocF):
            assocF = assoc
            F = Universe[i]
        if (assoc==0):
            Universe.pop(i)
    return [F,assocF,Universe]


def MMPC(TargetNode,Universe,X,alpha):
    CPC=[]
    print("Entering Phase I")
    print("MMPC_beggining: \n"+str(len(Universe)))
    while len(Universe)>0:
        CPC_old = list(CPC) # copy
        maxminheur=MaxMinHeuristic(TargetNode,CPC,list(Universe),X,alpha)
        F = maxminheur[0]
        assocF = maxminheur[1]
        Universe = maxminheur[2]
        if assocF > 0:
            CPC.append(F)
            indF=Universe.index(F)
            Universe.pop(indF)
        #if (len(CPC)==len(CPC_old)) or (len(CPC)>0.3*(len(Universe)-1)):
        if (len(CPC)==len(CPC_old)):
            break
        print("\nUniverse actual size:")
        print(len(Universe))
        print("CPC actual size:")
        print(len(CPC))
        print("CPC contents:")
        print_CPC(CPC)
    
    # Phase 2: Backward
    print("\nEntering Phase II")
    CPC=CPC[::-1]
    if len(CPC)>1:
        Z=list(CPC)
        for i in range(len(CPC)-1,-1,-1):
            # index is i
            Z.pop(i)
            #print("Analyzing D-separator for ",CPC[i]['name'])
            if ExistDseparator(TargetNode,CPC[i],Z,X,alpha)[0] == True:
                #print("it did exist! removing from cpc.")
                CPC.pop(i)
    return CPC

def arrayUniverse(TargetNode,arrayX):
    Universe = list(arrayX)
    for i in range(0,len(arrayX)):
        if (Universe[i]['name']==TargetNode):
            Universe.pop(i)
            break
    return Universe

def arrayX(X):
    returnable=[]
    for key in X:
        append_dict={'name':key,'data':X[key].copy(deep=True).tolist()}
        returnable.append(append_dict)
    return returnable

def MinAssoc(TargetNode, Xi,Z, fixedCondVars, X, alpha):
    
    min_assoc=999
    if len(Z)==0:
        min_assoc = Dep(TargetNode, Xi, fixedCondVars, X, alpha)
        subsetZ_min_assoc = fixedCondVars        
    else:
        #print(2**len(Z)-1)
        for IDsubsetZ_dec in range(1,2**len(Z)):
            IDsubsetZ_bin = bin(IDsubsetZ_dec)
            subsetZ = getZsubset(IDsubsetZ_bin,Z)            
            subsetZ_assoc=Dep(TargetNode, Xi, fixedCondVars+subsetZ,X,alpha)
            #counter+=1
            #print(subsetZ_assoc[IDsubsetZ_dec])
            if subsetZ_assoc < min_assoc:
                min_assoc = subsetZ_assoc
                subsetZ_min_assoc = fixedCondVars + subsetZ
                if (min_assoc==0):
                    break
    return min_assoc


def MMPC(TargetNode, X, alpha):
    # The universe will be an array of DataFrame Columns
    Universe = arrayUniverse(TargetNode,X)
    print(len(Universe))
    X = arrayX(X)
    print(len(X))
    for column in X:
        if (column['name']==TargetNode):
            TargetNode = column
            break
    CPC = MMPC_(TargetNode,Universe,X,alpha)
    return CPC