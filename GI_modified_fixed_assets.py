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

Total_Time = 300

a_file = open("SevenAssets.pkl", "rb")
assets = pickle.load(a_file)

a_file = open("InputPriceAndRHS.pkl", "rb")
scenarios = pickle.load(a_file)
Frc = scenarios['Prices'][:70]
RHSc = scenarios['RHS'][:70]
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
M.setParam('TimeLimit',Total_Time)
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
M.setParam('TimeLimit',Total_Time)
M.optimize()
M_scSlack = M
m_scr_slack = abs(M.ObjVal - M_scratch.ObjBound)/M.ObjVal*100

# Definition of function for GI

def Define_MILP_model(assets, RHSc, Frc, kn):
    KEYS = list(assets.keys())
    STORAGE = []; Z = []; SMP_T = []; SMP_P = []; SUMTURB = []; SUMPUMP = []; SPILLAGE = []; UT = []; UP = []
    #M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0                                                 = OptiSingle( RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0, p_T0, p_P0, u_T0, u_P0, v_T0, v_P0, w_T0, w_P0  = OptiSingle1(RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])
    
    XP  = M.addMVar(T, 0, 10**8)
    XT  = M.addMVar(T, 0, 10**8)
    OBJ = OBJ - np.ones(T)*10000@XP  - np.ones(T)*10000@XT
    M.setObjective(OBJ, GRB.MAXIMIZE)
    STORAGE.append(Storage0); Z.append(z0); 
    SMP_T.append(smp_T0); SMP_P.append(smp_P0)
    SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
    SPILLAGE.append(Spill0)
    UT.append(u_T0)
    UP.append(u_P0)
    for j in range(1, kn):
       #M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0                                                  = OptiSingle( RHSc, assets, KEYS[j], Frc, 60, env, RHSc*0, 0, M, OBJ)
        M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0, p_T0, p_P0, u_T0, u_P0, v_T0, v_P0, w_T0, w_P0  = OptiSingle1(RHSc, assets, KEYS[j], Frc, 60, env, RHSc*0, 0, M, OBJ)
    
        STORAGE.append(Storage0); Z.append(z0); 
        SMP_T.append(smp_T0); SMP_P.append(smp_P0)
        SUMTURB.append(SumTurb0); SUMPUMP.append(SumPump0)
        SPILLAGE.append(Spill0)
        UT.append(u_T0)
        UP.append(u_P0)
    M.addConstr(sum(SMP_T) + XT >= positiv(RHSc))
    M.addConstr(sum(SMP_P) + XP >= negativ(RHSc))
    return M, STORAGE, Z, SMP_T, SMP_P, SUMTURB, SUMPUMP, SPILLAGE, UT, UP

# Implementation of GI

def Gradual_Increase(assets, RHSc, Frc, TimeLimits, SubPools, u_TA, u_PA):
    M=[[]]*len(TimeLimits); STORAGE=[[]]*len(TimeLimits); Z=[[]]*len(TimeLimits); SMP_T=[[]]*len(TimeLimits); SMP_P=[[]]*len(TimeLimits); SUMTURB=[[]]*len(TimeLimits); SUMPUMP=[[]]*len(TimeLimits); SPILLAGE=[[]]*len(TimeLimits); UT=[[]]*len(TimeLimits); UP=[[]]*len(TimeLimits)
    m = len(TimeLimits)
    ii = 0
    M0, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model(assets, RHSc, Frc, SubPools[ii])
    M[ii] = M0; Z[ii] = Z0; STORAGE[ii] = STORAGE0; SMP_T[ii] = SMP_T[ii]; SMP_P[ii] = SMP_P0
    SUMTURB[ii] = SUMTURB0; SUMPUMP[ii] = SUMPUMP0; SPILLAGE[ii] = SPILLAGE0; UT[ii] = UT0; UP[ii] = UP0
    M[ii].setParam('TimeLimit',TimeLimits[ii])
    M[ii].optimize() 
    for ii in range(1, m):
        M0, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model(assets, RHSc, Frc, SubPools[ii])
        M[ii] = M0; Z[ii] = Z0; STORAGE[ii] = STORAGE0; SMP_T[ii] = SMP_T[ii]; SMP_P[ii] = SMP_P0
        SUMTURB[ii] = SUMTURB0; SUMPUMP[ii] = SUMPUMP0; SPILLAGE[ii] = SPILLAGE0; UT[ii] = UT0; UP[ii] = UP0
        
        for j in range(SubPools[ii-1]):
            Z[ii][j].start = Z[ii-1][j].x
            m = len(UT[ii-1][j])
            for i in range(m):
                UT[ii][j][i].start = UT[ii-1][j][i].x
                UP[ii][j][i].start = UP[ii-1][j][i].x
                
        
        for j in range(SubPools[ii-1], SubPools[ii]):
            m = len(UT[ii][j])
            for i in range(m):
                UT[ii][j][i].start = u_TA[j][i]
                UP[ii][j][i].start = u_PA[j][i]    
        M[ii].setParam('TimeLimit',TimeLimits[ii])
        M[ii].optimize() 
    return M[-1]

#Calculation of independent schedules
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


TimeLimits = [30, 270]
SubPools   = [3,   7]
M1 = Gradual_Increase(assets, RHSc, Frc, TimeLimits, SubPools, u_TA, u_PA)

TimeLimits = [50, 250]
SubPools   = [4,   7]
M2 = Gradual_Increase(assets, RHSc, Frc, TimeLimits, SubPools, u_TA, u_PA)

TimeLimits = [10,  20, 30, 240]
SubPools   = [2,   4,   6,   7]
M3 = Gradual_Increase(assets, RHSc, Frc, TimeLimits, SubPools, u_TA, u_PA)

print(' MIPgap of BF         = '+STR(M_scratch.MipGap*100)+'%')
print(' MIPgap of BF_modif   = '+STR(M_scSlack.MipGap*100)+'%')
print(' MIPgap of GI1        = '+STR(M1.MipGap*100)+'%')
print(' MIPgap of GI2        = '+STR(M2.MipGap*100)+'%')
print(' MIPgap of GI3        = '+STR(M3.MipGap*100)+'%')
END = timer()
print('It took '+STR(END - START)+' seconds to run this program')