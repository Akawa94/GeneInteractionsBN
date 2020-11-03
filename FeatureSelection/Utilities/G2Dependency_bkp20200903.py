import numpy as np
import pandas as pd
import scipy as sp
from scipy.stats import chi2
from math import floor
from math import ceil
import operator
from functools import reduce
import itertools
import random
import os
import time

class G2Dependency():
    
    def __init__(self, alpha, data_c_size):
        self._dep_dict = {}
        self._alpha = alpha
        self._data_c_size = data_c_size
        print('Succesful instantiation of Dep object')
        return
    
    def __get_cond_values(self, z, cond_vars):
        if len(cond_vars)>1:
            maxz = reduce(operator.mul,[int(x['name'].split('@')[1]) for x in cond_vars],1)
            d = np.zeros(len(cond_vars))

            div_dim = maxz/int(cond_vars[0]['name'].split('@')[1])
            num2div = z
            for i in range(0,len(cond_vars)-1):
                d[i] = floor(num2div/div_dim)
                num2div=num2div%div_dim
                div_dim = div_dim/int(cond_vars[i+1]['name'].split('@')[1])
            d[i+1] = num2div
        else:
            d = [z]
        return d 
    
    def __dep_dict_retriever(self, x_node, v_look):
        # making hash key
        x_node_key = x_node['name']+'_'+str(v_look)
        if x_node_key in self._dep_dict:
            return self._dep_dict[x_node_key]
        else:
            returnable = np.array(list(map(lambda x: 1 if x==v_look else 0, x_node['data'])))
            self._dep_dict[x_node_key] = returnable
            return returnable
    
    
    def put_arr(self,f_params,TargetNode,vvTarget,vvXi,Xi,CondVars):
        i,j,k = f_params
        flagData_arr = []
        condValue = self.__get_cond_values(k,CondVars)
        flagData_arr.append(np.ones(len(Xi['data'])))
        for l in range(0,len(CondVars)):
            X_l = CondVars[l]
            op_x = self.__dep_dict_retriever(X_l,condValue[l])
            flagData_arr.append(op_x)
        
        flagData_arr.append(self.__dep_dict_retriever(TargetNode,vvTarget[i]))
        flagData_arr.append(self.__dep_dict_retriever(Xi,vvXi[j]))
        
        return flagData_arr
    
    
    def dependency(self, TargetNode, Xi, CondVars, alpha):
        vvTarget = []
        for i in range(int(TargetNode['name'].split('@')[1])):
            vvTarget.append(i)

        vvXi = []
        for i in range(int(Xi['name'].split('@')[1])):
            vvXi.append(i)

        szCondVars = reduce(operator.mul,[int(x['name'].split('@')[1]) for x in CondVars],1)

        if self._data_c_size <= (5 * len(vvTarget)*len(vvXi)*(szCondVars)):
            return 1

        if (len(CondVars)==0): # test margianl dependency
            S = np.zeros((len(vvTarget),len(vvXi)))
            for i in range(0,len(vvTarget)):
                for j in range(0,len(vvXi)):
                    op1 = self.__dep_dict_retriever(TargetNode,vvTarget[i])
                    op2 = self.__dep_dict_retriever(Xi,vvXi[j])
                    S[i][j]=np.sum(op1*op2)
            G2=0
            N = np.sum(S)
            Si =np.sum(S,axis=1)
            Sj =np.sum(S,axis=0)
            Df = ((len(vvTarget)-1)*(len(vvXi)-1))
            Dedf = len((S[np.where(S>0)]))\
                          - len(Si[np.where(Si>0)])\
                          - len(Sj[np.where(Sj>0)]) +1
            if (Dedf<1):
                Dedf=1            
            for i in range(0,len(vvTarget)):
                for j in range(0,len(vvXi)):
                    if (S[i][j]>0):
                        G2 = G2 + S[i][j]*np.log((S[i][j])*N/(Si[i]*Sj[j]))

        else: # test conditional dependency
#             S = np.zeros((len(vvTarget),len(vvXi),szCondVars))
#             for i in range(0,len(vvTarget)):
#                 for j in range(0,len(vvXi)):
#                     for k in range(0,szCondVars):
#                         condValue = self.__get_cond_values(k,CondVars)
#                         flagData_arr = []
#                         flagData_arr.append(np.ones(len(Xi['data'])))

#                         for l in range(0,len(CondVars)):
#                             X_l = CondVars[l]

#                             op_x = self.__dep_dict_retriever(X_l,condValue[l])

#                             flagData_arr.append(op_x)
#                             #flagDataCondVars = flagDataCondVars*op_x

#                         #op1 = dep_dict_retriever(TargetNode,vvTarget[i])
#                         flagData_arr.append(self.__dep_dict_retriever(TargetNode,vvTarget[i]))
#                         #op2 = dep_dict_retriever(Xi,vvXi[j])
#                         flagData_arr.append(self.__dep_dict_retriever(Xi,vvXi[j]))

#                         flagDataCondVars = reduce(lambda x,y:x*y, flagData_arr)

#                         S[i][j][k]=np.sum(flagDataCondVars)
            i_t = len(vvTarget)
            j_t = len(vvXi)
            k_t = szCondVars
            mul_ijk = list(range(0,i_t*j_t*k_t))
            flagDataAcum = [self.put_arr([index//(j_t*k_t),index//(k_t)%j_t,index%k_t],TargetNode,vvTarget,vvXi,Xi,CondVars) for index in mul_ijk]
            flag_calc = *map(lambda x:np.sum(reduce(lambda y,z:y*z, np.array(x)),axis=0),flagDataAcum),
            S = np.asarray(flag_calc).reshape((i_t,j_t,k_t))
            
            
            G2 = 0
            Sjk = np.sum(S,axis=0)
            Sik = np.sum(S,axis=1)
            Sk = np.sum(Sjk,axis=0)

            Dedf = len(S[np.where(S>0)]) -\
                    len(Sik[np.where(Sik>0)]) -\
                    len(Sjk[np.where(Sjk>0)]) +\
                    len(Sk[np.where(Sk>0)])
            if Dedf<1:
                Dedf=1
                
            for i in range(0,len(vvTarget)):
                for j in range(0,len(vvXi)):
                    for k in range(0,szCondVars):
                        if (S[i][j][k]>0):
                            G2 = G2 + (S[i][j][k] * np.log(\
                                        (S[i][j][k] * Sk[k])/\
                                        (Sik[i][k] * Sjk[j][k])))
            G2 = 2*G2

        assoc = (alpha - (1 - chi2.cdf(G2,Dedf)))/alpha
        
        if assoc<0:
            assoc=0
            
        return assoc
        
             
