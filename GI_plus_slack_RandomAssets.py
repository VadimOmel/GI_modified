# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 22:44:42 2022

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
def subDict(flp, a1, a2):
    return {k: flp[k] for k in (a1, a2) if k in flp}

def Index_Of_Hydro_Asset_PP(Nn):
    return 'Hydro_Asset'+str(Nn)
def RandomNumber_Of_Turbines(m1,d):
    A=np.random.randint(m1, m1+d, 1)
    return A[0]

def Random_Int_Number(m1,d):
    A=np.random.randint(m1, m1+d, 1)
    return A[0]
    
def Random_Cont_Number(m1,d):
    A=np.random.uniform(m1, m1+d, 1)
    return A[0]

def GenerateTurbinesHydro(N_turbs):
    d3 = {}
    for j in range(N_turbs):
        Pmin = Random_Int_Number(200,600)
        Pmax = 2*Pmin
        Xbar = Random_Int_Number(50000,25000)
        cjT  = 2*Pmin*0.05
        bjT  = 2*Pmin*0.07
        d4 = {'Pmin': Pmin,
              'Pmax': 2*Pmin,
              'Xbar': Xbar,
              'cjT':  cjT,
              'bjT':  bjT,
              'ajT':  (Pmax - cjT - bjT)/Xbar,
              'StartCost': Random_Int_Number(2,11),
              'StopCost':  Random_Int_Number(2,11),
              'TminOn':    Random_Int_Number(2,6),
              'TminOff':   Random_Int_Number(2,7)}
        d3.update({'Turbine'+str(j+1) : d4})
    return d3

def Build_Hydro_Assets_on_Top(assets, Number_Additional_Assets):
    n1 = len(assets)
    for n in range(Number_Additional_Assets):
        N_turbs = Random_Int_Number(2,15)         
        assets.update({Index_Of_Hydro_Asset_PP(n+n1+1): 
                {
                    'Name': 'Asset'+str(n+n1+1),
                    'TechnicalFeature':
                    {
                    'Capacity': 2*Random_Cont_Number(2000000,4000000),    # kW
                    'N_Turbines': N_turbs,
                    'Turbines': GenerateTurbinesHydro(N_turbs)
        
                    }
                }
            })    
    return assets









a_file = open("SevenAssets.pkl", "rb")
assets = Build_Hydro_Assets_on_Top(assets, 12) # assets = pickle.load(a_file)

a_file = open("InputPriceAndRHS.pkl", "rb")
scenarios = pickle.load(a_file)
Frc = scenarios['Prices']
RHSc = scenarios['RHS']
T = len(Frc)

#SCRATCH
KEYS = list(assets.keys())
STORAGE = []; Z = []; SMP_T = []; SMP_P = []; SUMTURB = []; SUMPUMP = []; SPILLAGE = []
M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])

STORAGE.append(Storage0); Z.append(z0); 
SMP_T.append(smp_T0); SMP_P.append(smp_P0)
SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
SPILLAGE.append(Spill0)
for j in range(1, len(KEYS)):
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[j], Frc, 60, env, RHSc*0, 0, M, OBJ)
    STORAGE.append(Storage0); Z.append(z0); 
    SMP_T.append(smp_T0); SMP_P.append(smp_P0)
    SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
    SPILLAGE.append(Spill0)

M.addConstr(sum(SMP_T) >= positiv(RHSc))
M.addConstr(sum(SMP_P) >= negativ(RHSc))
M.setParam('TimeLimit',300)
M.optimize()
M_scratch = M

# SCRATCH CAP




KEYS = list(assets.keys())
STORAGE = []; Z = []; SMP_T = []; SMP_P = []; SUMTURB = []; SUMPUMP = []; SPILLAGE = []
M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])

XP  = M.addMVar(T, 0, 10**8)
XT  = M.addMVar(T, 0, 10**8)
OBJ = OBJ - np.ones(T)*10000@XP  - np.ones(T)*10000@XT

STORAGE.append(Storage0); Z.append(z0); 
SMP_T.append(smp_T0); SMP_P.append(smp_P0)
SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
SPILLAGE.append(Spill0)
for j in range(1, len(KEYS)):
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[j], Frc, 60, env, RHSc*0, 0, M, OBJ)
    STORAGE.append(Storage0); Z.append(z0); 
    SMP_T.append(smp_T0); SMP_P.append(smp_P0)
    SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
    SPILLAGE.append(Spill0)

M.addConstr(sum(SMP_T) + XT >= positiv(RHSc))
M.addConstr(sum(SMP_P) + XP >= negativ(RHSc))
M.setParam('TimeLimit',300)
M.optimize()
M_scSlack = M
m_scr_slack = abs(M.ObjVal - M_scratch.ObjBound)/M.ObjVal*100

kn = 3
KEYS = list(assets.keys())
STORAGE = []; Z = []; SMP_T = []; SMP_P = []; SUMTURB = []; SUMPUMP = []; SPILLAGE = []
M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])

XP  = M.addMVar(T, 0, 10**8)
XT  = M.addMVar(T, 0, 10**8)
OBJ = OBJ - np.ones(T)*10000@XP  - np.ones(T)*10000@XT
M.setObjective(OBJ, GRB.MAXIMIZE)
STORAGE.append(Storage0); Z.append(z0); 
SMP_T.append(smp_T0); SMP_P.append(smp_P0)
SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
SPILLAGE.append(Spill0)
for j in range(1, kn):
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[j], Frc, 60, env, RHSc*0, 0, M, OBJ)
    STORAGE.append(Storage0); Z.append(z0); 
    SMP_T.append(smp_T0); SMP_P.append(smp_P0)
    SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
    SPILLAGE.append(Spill0)

M.addConstr(sum(SMP_T) + XT >= positiv(RHSc))
M.addConstr(sum(SMP_P) + XP >= negativ(RHSc))

M.setParam('TimeLimit',85)
M.optimize()
Rt = M.RunTime
M.remove(M.getConstrs()[-2*T:])
for j in range(kn, len(KEYS)):
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[j], Frc, 60, env, RHSc*0, 0, M, OBJ)
    STORAGE.append(Storage0); Z.append(z0); 
    SMP_T.append(smp_T0); SMP_P.append(smp_P0)
    SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
    SPILLAGE.append(Spill0)
M.addConstr(sum(SMP_T) + XT >= positiv(RHSc))
M.addConstr(sum(SMP_P) + XP >= negativ(RHSc))
M.setParam('TimeLimit',300 - Rt)
M.optimize()

m4GI_gap = abs(M.ObjVal - M_scratch.ObjBound)/M.ObjVal*100

M_GI_4toAll = M

smpp = 0*SMP_T[0].x
smpp_p = 0*SMP_T[0].x
for j in range(len(SMP_T)):
    plt.plot(SMP_T[j].x)
    smpp = smpp + SMP_T[j].x
    smpp_p = smpp_p + SMP_P[j].x
j = 0
Visualize(SMP_T[j].x, Frc, 0, 200, 'prod', 'price', 'total production against price', 0)
j = j+1

Visualize(smpp, Frc, 0, 200, 'prod', 'price', 'total production against price', 0)

plt.plot(smpp)
plt.plot(positiv(RHSc))
plt.legend(('Cumulative Power','Minimum Required Power'), loc='upper left', shadow=True)
plt.ylim([0,75000])
plt.ylabel('Power (discharge/release)')
plt.xlabel('Time (h)')
plt.title('The case where classic GI cannot solve the problem')
plt.show()


plt.plot(smpp_p)
plt.plot(negativ(RHSc))
plt.legend(('Cumulative Power','Minimum Required Power'), loc='upper left', shadow=True)
plt.ylabel('Power (charge/pump)')
plt.xlabel('Time (h)')
plt.title('The case where classic GI cannot solve the problem')
plt.ylim([0,75000])
plt.show()



# SHUFFLING

from random import shuffle

KEYS = list(assets.keys())
shuffle(KEYS)

kn =2
STORAGE = []; Z = []; SMP_T = []; SMP_P = []; SUMTURB = []; SUMPUMP = []; SPILLAGE = []
M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])

XP  = M.addMVar(T, 0, 10**8)
XT  = M.addMVar(T, 0, 10**8)
OBJ = OBJ - np.ones(T)*10000@XP  - np.ones(T)*10000@XT
M.setObjective(OBJ, GRB.MAXIMIZE)
STORAGE.append(Storage0); Z.append(z0); 
SMP_T.append(smp_T0); SMP_P.append(smp_P0)
SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
SPILLAGE.append(Spill0)
for j in range(1, kn):
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[j], Frc, 60, env, RHSc*0, 0, M, OBJ)
    STORAGE.append(Storage0); Z.append(z0); 
    SMP_T.append(smp_T0); SMP_P.append(smp_P0)
    SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
    SPILLAGE.append(Spill0)

M.addConstr(sum(SMP_T) + XT >= positiv(RHSc))
M.addConstr(sum(SMP_P) + XP >= negativ(RHSc))

M.setParam('TimeLimit',20)
M.optimize()
Rt = M.RunTime
M.remove(M.getConstrs()[-2*T:])
for j in range(kn, len(KEYS)):
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[j], Frc, 60, env, RHSc*0, 0, M, OBJ)
    STORAGE.append(Storage0); Z.append(z0); 
    SMP_T.append(smp_T0); SMP_P.append(smp_P0)
    SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
    SPILLAGE.append(Spill0)
M.addConstr(sum(SMP_T) + XT >= positiv(RHSc))
M.addConstr(sum(SMP_P) + XP >= negativ(RHSc))
M.setParam('TimeLimit',300 - Rt)
M.optimize()
 
m4GI_gap_shuflfe = abs(M.ObjVal - M_scratch.ObjBound)/M.ObjVal*100
M_GI_4toAll_shuffle = M

print('MIPgap(       scratch       ) = '+STR(100*M_scratch.MIPgap)+'%')
print('MIPgap(    scratch_SLACK    ) = '+STR(100*M_scSlack.MIPgap)+'%')
print('MIPgap( GI [1-4]-All        ) = '+STR(100*M_GI_4toAll.MIPgap)+'%')
print('MIPgap( GI [1-4]-All shuffle) = '+STR(100*M_GI_4toAll_shuffle.MIPgap)+'%')
print('MIP gap with respect to the best upper bound')
print('MIPgap(       scratch       ) = '+STR(100*M_scratch.MIPgap)+'%')
print('MIPgap(    scratch_SLACK    ) = '+STR(m_scr_slack)+'%')
print('MIPgap( GI [1-4]-All        ) = '+STR(m4GI_gap)+'%')
print('MIPgap( GI [1-4]-All shuffle) = '+STR(m4GI_gap_shuflfe)+'%')


END = timer()
print('It took '+STR(END-START)+' seconds.')

kn = 13

StorageA = [[]]*len(assets); zA = [[]]*len(assets); smp_TA = [[]]*len(assets); smp_PA = [[]]*len(assets); SumTurbA = [[]]*len(assets); SumPumpA = [[]]*len(assets)
SpillA = [[]]*len(assets); p_TA = [[]]*len(assets); p_PA = [[]]*len(assets); u_TA = [[]]*len(assets); u_PA = [[]]*len(assets); v_TA = [[]]*len(assets); v_PA = [[]]*len(assets); w_TA = [[]]*len(assets); w_PA = [[]]*len(assets)

KEYS = list(assets.keys())
for j in range(0, len(KEYS)):
    M, OBJ, Storage, z, smp_T, smp_P, SumTurb, SumPump, Spill, p_T, p_P, u_T, u_P, v_T, v_P, w_T, w_P = OptiSingle1(RHSc*0.02, assets, KEYS[j], Frc, 10, env, RHSc*0, 1, [], [])
    StorageA[j] = Storage.x; zA[j] = z.x; smp_TA[j] = smp_T.x; smp_PA[j] = smp_P.x; SumTurbA[j]=SumTurb.x; SumPumpA[j]=SumPump.x; SpillA[j] = Spill.x
    m = len(p_T)
    for i in range(m):
        p_TA[j] = transform(p_T)
        u_TA[j] = transform(u_T)
        v_TA[j] = transform(v_T)
        w_TA[j] = transform(w_T)
        p_PA[j] = transform(p_P)
        u_PA[j] = transform(u_P)
        v_PA[j] = transform(v_P)
        w_PA[j] = transform(w_P)

KEYS = list(assets.keys())
STORAGE = []; Z = []; SMP_T = []; SMP_P = []; SUMTURB = []; SUMPUMP = []; SPILLAGE = []
M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])

XP  = M.addMVar(T, 0, 10**8)
XT  = M.addMVar(T, 0, 10**8)
OBJ = OBJ - np.ones(T)*10000@XP  - np.ones(T)*10000@XT
M.setObjective(OBJ, GRB.MAXIMIZE)
STORAGE.append(Storage0); Z.append(z0); 
SMP_T.append(smp_T0); SMP_P.append(smp_P0)
SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
SPILLAGE.append(Spill0)
for j in range(1, kn):
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0 = OptiSingle(RHSc, assets, KEYS[j], Frc, 60, env, RHSc*0, 0, M, OBJ)
    STORAGE.append(Storage0); Z.append(z0); 
    SMP_T.append(smp_T0); SMP_P.append(smp_P0)
    SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
    SPILLAGE.append(Spill0)

M.addConstr(sum(SMP_T) + XT >= positiv(RHSc))
M.addConstr(sum(SMP_P) + XP >= negativ(RHSc))

M.setParam('TimeLimit',85)
M.optimize()
Rt = M.RunTime
M.remove(M.getConstrs()[-2*T:])
for j in range(kn, len(KEYS)):
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0, p_T0, p_P0, u_T0, u_P0, v_T0, v_P0, w_T0, w_P0 = OptiSingle1(RHSc, assets, KEYS[j], Frc, 60, env, RHSc*0, 0, M, OBJ)
    STORAGE.append(Storage0); Z.append(z0); 
    SMP_T.append(smp_T0); SMP_P.append(smp_P0)
    SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
    SPILLAGE.append(Spill0)
    
    m = len(u_PA[j])
    for i in range(m):
        u_T0[i].start = u_TA[j][i]
        v_T0[i].start = v_TA[j][i]
        w_T0[i].start = w_TA[j][i]

        u_P0[i].start = u_PA[j][i]
        v_P0[i].start = v_PA[j][i]
        w_P0[i].start = w_PA[j][i]
    
M.addConstr(sum(SMP_T) + XT >= positiv(RHSc))
M.addConstr(sum(SMP_P) + XP >= negativ(RHSc))
M.setParam('TimeLimit',Total_Time - Rt)
M.optimize()

m4GI_gap = abs(M.ObjVal - M_scratch.ObjBound)/M.ObjVal*100

M_GI_4toAll = M