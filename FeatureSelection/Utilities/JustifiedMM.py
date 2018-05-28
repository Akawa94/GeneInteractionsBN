def MaxMinHeuristic(TargetNode, CPC, Universe, X, alpha):
    F_arr=[]
    Z = CPC
    fixedCondVars = []
    if (len(CPC)>0):
        Z = CPC[0:-1]           # all but the last one   
        fixedCondVars = [CPC[-1]] # we use last one
    for i in range(len(Universe)-1,-1,-1):
        if (len(Universe[i])==0):
            continue
        assoc = MinAssoc(TargetNode,Universe[i],Z,fixedCondVars,X,alpha)
        if (assoc>0):
            F_arr.append([Universe[i],assoc])
        if (assoc==0):
            Universe.pop(i)
            
    return [F_arr,Universe]

def MMPC_(TargetNode,Universe,X,alpha):
    
    CPC=[]    
    # Phase I: Foward
    print("Entering Phase I")
    print("MMPC_beggining: \n"+str(len(Universe)))
    while len(Universe)>0:
        CPC_old = list(CPC) # copy
        maxminheur=MaxMinHeuristic(TargetNode,CPC,list(Universe),X,alpha)
        F_arr = maxminheur[0]
        if (len(F_arr)>0):
            max_assoc = max([x[1] for x in F_arr])
        else:
            max_assoc = 0
        
        Universe = maxminheur[1]
        if max_assoc > 0:
            # append all nodes that fulfill the requirements
            for nn in [x[0] for x in F_arr if x[1]>=max_assoc]:
                CPC.append(nn)
                indF=Universe.index(nn)
                Universe.pop(indF)
        
        if (len(CPC)==len(CPC_old)):
            break
        print("\nUniverse actual size:")
        print(len(Universe))
        print("CPC size before filtering: "+str(len(CPC)))
        # Embedded RW in backward_phase
        CPC = backward_phase_v4(TargetNode,CPC,X,alpha,len(CPC))
        print("CPC size after filtering: "+str(len(CPC)))
        print(list(map(lambda x:x['name'],CPC)))
        
    return CPC

def backward_phase_v2(TargetNode, CPC, X,alpha):
    #print("\nEntering Phase II")
    if (len(CPC)>1):
        # will divide the sets in 
        Z=list(CPC)
        while len(Z)>0:
            # Will test all vs the Target, and select the lowest
            minDep=999
            minDepIndex=-1
            for i in range(0,len(Z)):
                zCopy = list(Z)
                zCopy.pop(i)
                auxDep = Dep(TargetNode,Z[i],zCopy,X,alpha)
                if (minDep>auxDep):
                    minDep = auxDep
                    minDepIndex=i
            # Find index in CPC
            CPCindex=-1
            for i in range(0,len(CPC)):
                if (CPC[i]['name']==Z[minDepIndex]['name']):
                    CPCindex=i
                    break
            
            # Remove from Z
            Z.pop(minDepIndex)
            #print_names(Z)
            #print("Analyzing D-separator for ",CPC[CPCindex]['name'])
            if ExistDseparator(TargetNode,CPC[CPCindex],Z,X,alpha)[0] == True:
                #print("it did exist! removing from cpc.")
                CPC.pop(CPCindex)
    return CPCs


def backward_phase_v4(TargetNode, CPC, X,alpha,loops_flag):
    #print("\nEntering Phase II")
    votation_arr = []
    if (len(CPC)>1):
        # will divide the sets in np.log(size)*2
        size_sets = ceil(np.log(len(CPC))*2)
        #print(size_sets)
        n_sets = ceil(len(CPC)/size_sets)
        for k in range(0,loops_flag):
            nodes_sets = []
            CPC_copy = list(CPC)
            for set_counter in range(0,n_sets):
                sub_node_set = []
                while (len(sub_node_set)<size_sets and len(CPC_copy)>0):
                    rand_index = random.randint(0,len(CPC_copy))-1
                    sub_node_set.append(CPC_copy[rand_index])
                    CPC_copy.pop(rand_index)
                nodes_sets.append(sub_node_set)
            pre_cpc = []
            for n_set in nodes_sets:
                for node in backward_phase_v2(TargetNode,n_set,X,alpha):
                    pre_cpc.append(node)
            votation_arr.append(pre_cpc)
        # voting
        returnable_cpc = []
        for node in CPC:
            v_counter=0
            for votation in votation_arr:
                if node in votation:
                    v_counter=v_counter+1
            if (v_counter>(loops_flag-loops_flag**(1/2))):
                returnable_cpc.append(node)
        return returnable_cpc
    else:
        return CPC


def backward_phase_v2(TargetNode, CPC, X,alpha):
    #print("\nEntering Phase II")
    if (len(CPC)>1):
        # will divide the sets in 
        Z=list(CPC)
        while len(Z)>0:
            # Will test all vs the Target, and select the lowest
            minDep=999
            minDepIndex=-1
            for i in range(0,len(Z)):
                zCopy = list(Z)
                zCopy.pop(i)
                auxDep = Dep(TargetNode,Z[i],zCopy,X,alpha)
                if (minDep>auxDep):
                    minDep = auxDep
                    minDepIndex=i
            # Find index in CPC
            CPCindex=-1
            for i in range(0,len(CPC)):
                if (CPC[i]['name']==Z[minDepIndex]['name']):
                    CPCindex=i
                    break
            
            # Remove from Z
            Z.pop(minDepIndex)
            #print_names(Z)
            #print("Analyzing D-separator for ",CPC[CPCindex]['name'])
            if ExistDseparator(TargetNode,CPC[CPCindex],Z,X,alpha)[0] == True:
                #print("it did exist! removing from cpc.")
                CPC.pop(CPCindex)
    return CPC