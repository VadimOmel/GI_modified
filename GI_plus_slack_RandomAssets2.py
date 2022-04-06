# -*- coding: utf-8 -*-
"""
Created on Thu Feb 10 09:25:05 2022

@author: OMELCHENKO_V
"""

 
import copy
import pickle
import numpy as np
import gurobipy as gp
from timeit import default_timer as timer
import matplotlib.pyplot as plt
timeStart = timer()
from q_Opti3 import *
from gurobipy import GRB, Model, Env
START = timer()

connection_params = {
   "ComputeServer":  "YOUR_SERVER",
   "Username": "YOUR_USERNAME",
   "ServerPassword": "YOUR_PASSWORD"
    }

env = gp.Env(params=connection_params)






N_assets    = 50
Total_Time  = 500

a_file = open("SevenAssets.pkl", "rb")
assets = Build_Hydro_Assets_on_Top({}, N_assets) # assets = pickle.load(a_file)

a_file = open("InputPriceAndRHS1.pkl", "rb")
scenarios = pickle.load(a_file)
Frc = scenarios['Prices']
RHSc = scenarios['RHS']
T = len(Frc)

#scratch
M0, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model0(assets, RHSc, Frc, N_assets, env)
M0.setParam('OutputFlag', 1)
M0.setParam('TimeLimit',Total_Time)
M0.optimize()

#scratch + slack var
M1, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model(assets, RHSc, Frc, N_assets, env)
M1.setParam('OutputFlag', 1)
M1.setParam('TimeLimit',Total_Time)
M1.optimize()





# optimization of each PP separately. Can be conducted in parallel
StorageA = [[]]*len(assets); zA = [[]]*len(assets); smp_TA = [[]]*len(assets); smp_PA = [[]]*len(assets); SumTurbA = [[]]*len(assets); SumPumpA = [[]]*len(assets)
SpillA = [[]]*len(assets); p_TA = [[]]*len(assets); p_PA = [[]]*len(assets); u_TA = [[]]*len(assets); u_PA = [[]]*len(assets); v_TA = [[]]*len(assets); v_PA = [[]]*len(assets); w_TA = [[]]*len(assets); w_PA = [[]]*len(assets)

KEYS = list(assets.keys())
UpperBoundConstraintFree = 0
RUNTIMES = []
OBJVALS  = []
SMP_TT   = []
SMP_PP   = []
for j in range(0, len(KEYS)):
    
    M, OBJ, Storage, z, smp_T, smp_P, SumTurb, SumPump, Spill, p_T, p_P, u_T, u_P, v_T, v_P, w_T, w_P = OptiSingle1(RHSc*0, assets, KEYS[j], Frc, 6, env, RHSc*0, 1, [], [])
    StorageA[j] = Storage.x; zA[j] = z.x; smp_TA[j] = smp_T.x; smp_PA[j] = smp_P.x; SumTurbA[j]=SumTurb.x; SumPumpA[j]=SumPump.x; SpillA[j] = Spill.x
    m = len(p_T)
    SMP_TT.append(smp_T.x)
    SMP_PP.append(smp_P.x)
    RUNTIMES.append(M.RunTime)
    OBJVALS.append(M.ObjVal)
    UpperBoundConstraintFree += M.ObjVal
    print('asset = '+KEYS[j]+', last runTime = '+STR(M.RunTime)+', current ub = '+STR(UpperBoundConstraintFree))
    for i in range(m):
        p_TA[j] = transform(p_T)
        u_TA[j] = transform(u_T)
        v_TA[j] = transform(v_T)
        w_TA[j] = transform(w_T)
        p_PA[j] = transform(p_P)
        u_PA[j] = transform(u_P)
        v_PA[j] = transform(v_P)
        w_PA[j] = transform(w_P)



#GI+slack
TL1 = 150
TimeLimits = [TL1, Total_Time - TL1]
SubPools   = [30,          N_assets]
M2, smp_T, smp_P, Z_storage, storage, ut, up = Gradual_Increase(assets, RHSc, Frc, TimeLimits, SubPools, u_TA, u_PA, env)

m_objval   =  max([M0.ObjVal,   M1.ObjVal,    M2.ObjVal  ])
m_objbound =  min([M0.ObjBound, M1.ObjBound,  M2.ObjBound])

mipgap = abs(m_objval-m_objbound)/m_objval

print('Writing out results: MIPgap, ObjValues, ObjBounds')
print('Best MIPgap:    scr = '+STR(100*M0.MIPgap)  +'%,       scr+slack = '+STR(100*M1.MIPgap)  +'%,       GI+slack = ' +STR(100*M2.MIPgap)  +'%,       best val = ' +STR(100*mipgap )+'%'  )
print('Best ObjBound:  scr = '+STR(M0.ObjBound)+', scr+slack = '+STR(M1.ObjBound)+', GI+slack = ' +STR(M2.ObjBound)+', best val = ' +STR(m_objbound))
print('Best Objective: scr = '+STR(M0.ObjVal)  +', scr+slack = '+STR(M1.ObjVal)  +', GI+slack = ' +STR(M2.ObjVal)  +', best val = ' +STR(m_objval)  )
print('Maximum possible value: '+STR(UpperBoundConstraintFree))


# abs(M.ObjVal - M_scratch.ObjBound)/M.ObjVal*100

def CumulativeProduction0(SMP_T, SMP_P):
    smpp    =  0*SMP_T[0]
    smpp_p  =  0*SMP_T[0]
    for j in range(len(SMP_T)):
        #plt.plot(SMP_T[j].x)
        smpp   = smpp   + SMP_T[j]
        smpp_p = smpp_p + SMP_P[j]    
    return smpp, smpp_p 
    
def CumulativeProduction(SMP_T, SMP_P):
    smpp    =  0*SMP_T[0].x
    smpp_p  =  0*SMP_T[0].x
    for j in range(len(SMP_T)):
        #plt.plot(SMP_T[j].x)
        smpp   = smpp   + SMP_T[j].x
        smpp_p = smpp_p + SMP_P[j].x    
    return smpp, smpp_p

smpp,   smpp_p = CumulativeProduction( smp_T,  smp_P )
smpp1, smpp_p1 = CumulativeProduction0(SMP_TT, SMP_PP)

Visualize(smpp, Frc, 0, 200, 'prod', 'price', 'total production (release) against price', 0)

Visualize(smpp_p, -Frc, 0, 200, 'prod', 'price*(-1)', 'total production (pump) against price*(-1)', 0)

plt.plot(smpp1)
plt.plot(smpp)
plt.legend(('TotPower (No coupling)','TotPower (Coupling)'), loc='upper left', shadow=True)
plt.title('Total Power (Released)')
plt.show()

plt.plot(smpp_p1)
plt.plot(smpp_p)
plt.ylim([0,1.3*max(smpp_p1)])
plt.legend(('TotPower (No coupling)','TotPower (Coupling)'), loc='upper left', shadow=True)
plt.title('Total Power (Pumped)')
plt.show()

plt.plot(smpp1)
plt.plot(smpp)
plt.plot(RHSc)
plt.legend(('TotPower (No coupling)','TotPower (Coupling)','Task'), loc='upper left', shadow=True)
plt.ylim([0,1.34*max(smpp)])
plt.title('Total Power (Pumped)')
plt.show()

plt.plot(OBJVALS)
plt.title('The optimal value of each asset - independently')
plt.show()

plt.plot(smpp)
plt.plot(positiv(RHSc))
plt.legend(('Cumulative Power','Minimum Required Power'), loc='upper left', shadow=True)
plt.ylim([0,1.3*max(smpp)])
plt.ylabel('Power (discharge/release)')
plt.xlabel('Time (h)')
plt.title('The case where classic GI cannot solve the problem')
plt.show()


# plt.plot(smpp_p)
# plt.plot(negativ(RHSc))
# plt.legend(('Cumulative Power','Minimum Required Power'), loc='upper right', shadow=True)
# plt.ylabel('Power (charge/pump)')
# plt.xlabel('Time (h)')
# plt.title('The case where classic GI cannot solve the problem')
# #plt.ylim([0,75000])
# plt.show()

#scratch
# M00, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model0(assets, RHSc*0, Frc, 5, env)
# M00.setParam('OutputFlag', 1)
# M00.setParam('TimeLimit',30)
# M00.optimize()
# #GI (classic)
# TimeLimits = [30, 270]
# SubPools   = [18,   N_assets]
# M3 = Gradual_Increase0(assets, RHSc, Frc, TimeLimits, SubPools, u_TA, u_PA, env)
# # Scratch + slack
# M11, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model(assets, RHSc, Frc, N_assets, env)
# M11.setParam('TimeLimit',300)
# M11.optimize()
# # Scratch
# M01, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model0(assets, RHSc, Frc, N_assets, env)
# M01.setParam('TimeLimit',300)
# M01.optimize()
# #GI+slack
# TL1 = 50
# TimeLimits = [TL1, 300 - TL1]
# SubPools   = [12,          N_assets]
# M02 = Gradual_Increase(assets, RHSc, Frc, TimeLimits, SubPools, u_TA, u_PA, env)



# m_objval1   =  max([M01.ObjVal,   M11.ObjVal,    M02.ObjVal,    M3.ObjVal  ])
# m_objbound1 =  min([M01.ObjBound, M11.ObjBound,  M02.ObjBound, M3.ObjBound])

# mipgap1 = abs(m_objval1-m_objbound1)/m_objval1

# print([M3.MIPgap, M11.MIPgap, M01.MIPgap, M02.MIPgap, mipgap1])
