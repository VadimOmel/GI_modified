# -*- coding: utf-8 -*-
"""
Created on Wed Dec 29 13:17:12 2021

@author: OMELCHENKO_V
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 18:32:00 2021

@author: OMELCHENKO_V
"""

# -*- coding: utf-8 -*-
"""
Created on Tue Dec 28 00:53:53 2021

@author: OMELCHENKO_V
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 12:23:51 2021

@author: OMELCHENKO_V
"""

# -*- coding: utf-8 -*-
"""
Created on Sat Dec 25 18:59:23 2021

@author: OMELCHENKO_V
"""

from timeit import default_timer as timer
import numpy as np
import copy
import matplotlib.pyplot as plt
timeStart = timer()
import gurobipy as gp
from gurobipy import GRB, Model 

def positiv(vector):
    vect = []
    for v in vector:
        if v >= 0:
            vect.append(v)
        else:
            vect.append(0)
    return np.array(vect)

def negativ(vector):
    vect = []
    for v in vector:
        if v <= 0:
            vect.append(-v)
        else:
            vect.append(0)
    return np.array(vect)



def CUMSUMV(vector, N):
    n = len(vector)
    if N<=n:
        N = n
    VcOutp = []
    for j in range(n):
        VcOutp.append(sum(vector[j:]))
    for j in range(n+1,N+1):
        VcOutp.append(0)
    return VcOutp
def OptiSubsetGIturbine(RHSc, assets, Frc1, Frc2, TL, env):
    Frc = Frc1
    T = len(Frc)
    L0 = 0
    PowerPl = 'HydroPowerPlant1'
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    Penlt = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['StartCost']
    
    Inflows = np.random.randint(2, size=len(Frc))
    Inflows = Inflows*Capacity/(180*sum(Inflows))
    bStorage = np.cumsum(Inflows*0)
    bStorage = L0 + np.array(bStorage)
    
    M = Model(env=env)
    Storage  = M.addMVar(T, 0, Capacity)
    Spill  = M.addMVar(T, 0, Capacity)
    z      = M.addMVar(T, 0, 1)
    z.vType = gp.GRB.BINARY
    
    SumPmax  = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    SmPowerMax1 = SumPmax1
    SumTurb = M.addMVar(T, 0, SumPmax)
    SumPump = M.addMVar(T, 0, SumPmax)
    
    smp_T   = M.addMVar(T, 0, SumPmax1)
    smp_P   = M.addMVar(T, 0, SumPmax1)
    
    
    
    
    q        = []
    p_T      = []
    u_T      = []
    v_T      = []
    w_T      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        q.append(M.addMVar(T, 0, assets['HydroPowerPlant1']['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T.append(M.addMVar(T, 0, 1))
        v_T.append(M.addMVar(T, 0, 1))
        w_T.append(M.addMVar(T, 0, 1))
        p_T.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumTurb - sum(q)   == 0)
    M.addConstr(smp_T   - sum(p_T) == 0)
    
    
    qpump = []
    p_P     = []
    u_P      = []
    v_P      = []
    w_P      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump.append(M.addMVar(T, 0, assets['HydroPowerPlant1']['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P.append(M.addMVar(T, 0, 1))
        v_P.append(M.addMVar(T, 0, 1))
        w_P.append(M.addMVar(T, 0, 1))
        p_P.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump - sum(qpump) == 0)
    M.addConstr(smp_P   - sum(p_P)   == 0)
    
    UT       = assets['HydroPowerPlant1']['TechnicalFeature']['Turbines']['Turbine1']['TminOn:']
    DT       = assets['HydroPowerPlant1']['TechnicalFeature']['Turbines']['Turbine1']['TminOff:']
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    M.addConstr(IDD@Storage + A1@SumTurb - A1@SumPump + A1@Spill == bStorage)
    UU = np.zeros(len(p_T))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T)):
        M.addConstr(IDD@u_T[jj] - A1@v_T[jj] + A1@w_T[jj] == UU[jj])
        M.addConstr(IDD@u_P[jj] - A1@v_P[jj] + A1@w_P[jj] == UU[jj])
        M.addConstr(IDD@u_T[jj]+IDD@u_P[jj] <= 1)
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P[jj]-aj*IDD@qpump[jj]          >= 0 )
        M.addConstr(IDD@p_T[jj]-aj*IDD@q[jj] - cj*IDD @ z <= bj)
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
        
        M.addConstr(IDD @ p_P[jj] - Pmin01*IDD @ u_P[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_P[jj] - IDD @ p_P[jj] >= 0)
        
        M.addConstr(IDD @ p_T[jj] - Pmin01*IDD @ u_T[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_T[jj] - IDD @ p_T[jj] >= 0)
    
        M.addConstr(IDD@u_T[jj] - BU@v_T[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T[jj] + BT@w_T[jj] <=onesVec-wAdd)
    
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          >= 0 )
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          <= Capacity/2 )
        
    
    
    PowerPl = 'HydroPowerPlant2'
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    
    Inflows = np.random.randint(2, size=len(Frc))
    Inflows = Inflows*Capacity/(180*sum(Inflows))
    bStorage = np.cumsum(Inflows*0)
    bStorage = L0 + np.array(bStorage)
    
    
    Storage1  = M.addMVar(T, 0, Capacity)
    Spill1  = M.addMVar(T, 0, Capacity)
    z1      = M.addMVar(T, 0, 1)
    z1.vType = gp.GRB.BINARY
    
    SumPmax = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    SmPowerMax2 = SumPmax1
    SumTurb1 = M.addMVar(T, 0, SumPmax)
    SumPump1 = M.addMVar(T, 0, SumPmax)
    smp_T1   = M.addMVar(T, 0, SumPmax1)
    smp_P1   = M.addMVar(T, 0, SumPmax1)
    
    q1        = []
    p_T1      = []
    u_T1      = []
    v_T1      = []
    w_T1      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        q1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T1.append(M.addMVar(T, 0, 1))
        v_T1.append(M.addMVar(T, 0, 1))
        w_T1.append(M.addMVar(T, 0, 1))
        p_T1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T1[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumTurb1 - sum(q1)   == 0)
    M.addConstr(smp_T1   - sum(p_T1) == 0)
    
    qpump1 = []
    p_P1     = []
    u_P1      = []
    v_P1      = []
    w_P1      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P1.append(M.addMVar(T, 0, 1))
        v_P1.append(M.addMVar(T, 0, 1))
        w_P1.append(M.addMVar(T, 0, 1))
        p_P1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P1[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump1 - sum(qpump1) == 0)
    M.addConstr(smp_P1    - sum(p_P1) == 0)
    
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    M.addConstr(IDD@Storage1 + A1@SumTurb1 - A1@SumPump1 - A1@SumTurb + A1@SumPump + A1@Spill1 == bStorage)
    UU = np.zeros(len(p_T1))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T1)):
        M.addConstr(IDD@u_T1[jj] - A1@v_T1[jj] + A1@w_T1[jj] == UU[jj])
        M.addConstr(IDD@u_P1[jj] - A1@v_P1[jj] + A1@w_P1[jj] == UU[jj])
        M.addConstr(IDD@u_T1[jj]+IDD@u_P1[jj] <= 1)
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P1[jj]-aj*IDD@qpump1[jj]          >= 0 )
        M.addConstr(IDD@p_T1[jj]-aj*IDD@q1[jj] - cj*IDD @ z1 <= bj)
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
        
        
        M.addConstr(IDD @ p_P1[jj] - Pmin01*IDD @ u_P1[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_P1[jj] - IDD @ p_P1[jj] >= 0)
        
        M.addConstr(IDD @ p_T1[jj] - Pmin01*IDD @ u_T1[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_T1[jj] - IDD @ p_T1[jj] >= 0)
    
        M.addConstr(IDD@u_T1[jj] - BU@v_T1[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T1[jj] + BT@w_T1[jj] <=onesVec-wAdd)
    
    M.addConstr(IDD @ Storage1-Capacity/2 * IDD @ z1          >= 0 )
    M.addConstr(IDD @ Storage1-Capacity/2 * IDD @ z1          <= Capacity/2 )
        
    
    #RHSc = np.sin(np.cumsum(np.ones(T))/20)*4000
    M.addConstr(IDD @ smp_T + IDD @ smp_T1 >= RHSc)
    M.addConstr(IDD @ smp_T + IDD @ smp_T1 >= 0.1*negativ(RHSc))

        
    positID = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    negatID = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    
    
    M.addConstr(IDD @ positID - positiv(RHSc) == 0)
    M.addConstr(IDD @ negatID - 0.1*negativ(RHSc) == 0)
    
    positDA = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    negatDA = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    
    M.addConstr(IDD @ positDA - IDD @ smp_T - IDD @ smp_T1 + positiv(RHSc) == 0)
    M.addConstr(IDD @ negatDA - IDD @ smp_P - IDD @ smp_P1 + 0.1*negativ(RHSc) == 0)
    
    
        
    OBJ = 0
    OBJ = OBJ + Frc @ positDA - Frc @ negatDA
    for sc in range(len(p_T)):
        OBJ = OBJ - Penlt*np.ones(T) @ v_T[sc] - Penlt*np.ones(T) @ w_T[sc]\
            - Penlt*np.ones(T)@v_P[sc] - Penlt*np.ones(T)@w_P[sc] 
    
    
    for sc in range(len(p_T1)):
        OBJ = OBJ - Penlt*np.ones(T)@v_T1[sc] - Penlt*np.ones(T)@w_T1[sc]\
            - Penlt*np.ones(T)@v_P1[sc] - Penlt*np.ones(T)@w_P1[sc] 
    OBJ = OBJ - Penlt*np.ones(T) @ Spill
    OBJ = OBJ - Penlt*np.ones(T) @ Spill1
    OBJ = OBJ + Frc2 @ positID - Frc2 @ negatID
    
    M.setParam('TimeLimit',TL)
    M.setObjective(OBJ, GRB.MAXIMIZE)
    M.optimize()
    ppT = []
    uuT = []
    vvT = []
    wwT = []
    qq  = []
    for jj in range(len(p_T)):
        ppT.append(p_T[jj].x)
        uuT.append(u_T[jj].x)
        vvT.append(v_T[jj].x)
        wwT.append(w_T[jj].x)
        qq.append(q[jj].x)
    ppP = []
    uuP = []
    vvP = []
    wwP = []
    qpm  = []
    for jj in range(len(p_P)):
        ppP.append(p_P[jj].x)
        uuP.append(u_P[jj].x)
        vvP.append(v_P[jj].x)
        wwP.append(w_P[jj].x)
        qpm.append(qpump[jj].x)
    Strg = Storage.x
    Spll = Spill.x
    zzzz = z.x
    SmTrb = SumTurb.x
    SmPmp = SumPump.x
    
    ppT1 = []
    uuT1 = []
    vvT1 = []
    wwT1 = []
    qq1  = []
    for jj in range(len(p_T1)):
        ppT1.append(p_T1[jj].x)
        uuT1.append(u_T1[jj].x)
        vvT1.append(v_T1[jj].x)
        wwT1.append(w_T1[jj].x)
        qq1.append(q1[jj].x)
    ppP1 = []
    uuP1 = []
    vvP1 = []
    wwP1 = []
    qpm1  = []
    for jj in range(len(p_P1)):
        ppP1.append(p_P1[jj].x)
        uuP1.append(u_P1[jj].x)
        vvP1.append(v_P1[jj].x)
        wwP1.append(w_P1[jj].x)
        qpm1.append(qpump1[jj].x)
    Strg1 = Storage1.x
    Spll1 = Spill1.x
    zzzz1 = z1.x
    SmTrb1 = SumTurb1.x
    SmPmp1 = SumPump1.x
    return Strg,  Spll,  zzzz,  SmTrb,  SmPmp,\
           ppT,   uuT,   vvT,   wwT,    qq, \
           ppP,   uuP,   vvP,   wwP,    qpm,\
           Strg1, Spll1, zzzz1, SmTrb1, SmPmp1,\
           ppT1,  uuT1,  vvT1,  wwT1,   qq1, \
           ppP1,   uuP1,  vvP1,  wwP1,   qpm1, M
    




def OptiEntireWS(RHSc, assets, Frc1, Frc2, TL, env, Strg,  Spll,  zzzz,  SmTrb,  SmPmp, ppT,  uuT,   vvT,   wwT,    qq, ppP,   uuP,   vvP,   wwP,    qpm, Strg1, Spll1, zzzz1, SmTrb1, SmPmp1, ppT1,  uuT1,  vvT1,  wwT1,   qq1, ppP1,   uuP1,  vvP1,  wwP1,   qpm1):
    Frc = Frc1
    T = len(Frc)
    L0 = 0
    PowerPl = 'HydroPowerPlant1'
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    Penlt = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['StartCost']
    
    Inflows = np.random.randint(2, size=len(Frc))
    Inflows = Inflows*Capacity/(180*sum(Inflows))
    bStorage = np.cumsum(Inflows*0)
    bStorage = L0 + np.array(bStorage)
    
    M = Model(env=env)
    Storage  = M.addMVar(T, 0, Capacity)
    Spill  = M.addMVar(T, 0, Capacity)
    z      = M.addMVar(T, 0, 1)
    z.vType = gp.GRB.BINARY
    
    SumPmax  = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    SmPowerMax1 = SumPmax1
    SumTurb = M.addMVar(T, 0, SumPmax)
    SumPump = M.addMVar(T, 0, SumPmax)
    
    smp_T   = M.addMVar(T, 0, SumPmax1)
    smp_P   = M.addMVar(T, 0, SumPmax1)
    
    
    
    
    q        = []
    p_T      = []
    u_T      = []
    v_T      = []
    w_T      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        q.append(M.addMVar(T, 0, assets['HydroPowerPlant1']['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T.append(M.addMVar(T, 0, 1))
        v_T.append(M.addMVar(T, 0, 1))
        w_T.append(M.addMVar(T, 0, 1))
        p_T.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumTurb - sum(q)   == 0)
    M.addConstr(smp_T   - sum(p_T) == 0)
    
    
    qpump = []
    p_P     = []
    u_P      = []
    v_P      = []
    w_P      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump.append(M.addMVar(T, 0, assets['HydroPowerPlant1']['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P.append(M.addMVar(T, 0, 1))
        v_P.append(M.addMVar(T, 0, 1))
        w_P.append(M.addMVar(T, 0, 1))
        p_P.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump - sum(qpump) == 0)
    M.addConstr(smp_P   - sum(p_P)   == 0)
    
    UT       = assets['HydroPowerPlant1']['TechnicalFeature']['Turbines']['Turbine1']['TminOn:']
    DT       = assets['HydroPowerPlant1']['TechnicalFeature']['Turbines']['Turbine1']['TminOff:']
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    M.addConstr(IDD@Storage + A1@SumTurb - A1@SumPump + A1@Spill == bStorage)
    UU = np.zeros(len(p_T))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T)):
        M.addConstr(IDD@u_T[jj] - A1@v_T[jj] + A1@w_T[jj] == UU[jj])
        M.addConstr(IDD@u_P[jj] - A1@v_P[jj] + A1@w_P[jj] == UU[jj])
        M.addConstr(IDD@u_T[jj]+IDD@u_P[jj] <= 1)
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P[jj]-aj*IDD@qpump[jj]          >= 0 )
        M.addConstr(IDD@p_T[jj]-aj*IDD@q[jj] - cj*IDD @ z <= bj)
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
        
        M.addConstr(IDD @ p_P[jj] - Pmin01*IDD @ u_P[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_P[jj] - IDD @ p_P[jj] >= 0)
        
        M.addConstr(IDD @ p_T[jj] - Pmin01*IDD @ u_T[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_T[jj] - IDD @ p_T[jj] >= 0)
    
        M.addConstr(IDD@u_T[jj] - BU@v_T[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T[jj] + BT@w_T[jj] <=onesVec-wAdd)
    
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          >= 0 )
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          <= Capacity/2 )
        
    
    
    PowerPl = 'HydroPowerPlant2'
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    
    Inflows = np.random.randint(2, size=len(Frc))
    Inflows = Inflows*Capacity/(180*sum(Inflows))
    bStorage = np.cumsum(Inflows*0)
    bStorage = L0 + np.array(bStorage)
    
    
    Storage1  = M.addMVar(T, 0, Capacity)
    Spill1  = M.addMVar(T, 0, Capacity)
    z1      = M.addMVar(T, 0, 1)
    z1.vType = gp.GRB.BINARY
    
    SumPmax = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    SmPowerMax2 = SumPmax1
    SumTurb1 = M.addMVar(T, 0, SumPmax)
    SumPump1 = M.addMVar(T, 0, SumPmax)
    smp_T1   = M.addMVar(T, 0, SumPmax1)
    smp_P1   = M.addMVar(T, 0, SumPmax1)
    
    q1        = []
    p_T1      = []
    u_T1      = []
    v_T1      = []
    w_T1      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        q1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T1.append(M.addMVar(T, 0, 1))
        v_T1.append(M.addMVar(T, 0, 1))
        w_T1.append(M.addMVar(T, 0, 1))
        p_T1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T1[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumTurb1 - sum(q1)   == 0)
    M.addConstr(smp_T1   - sum(p_T1) == 0)
    
    qpump1 = []
    p_P1     = []
    u_P1      = []
    v_P1      = []
    w_P1      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P1.append(M.addMVar(T, 0, 1))
        v_P1.append(M.addMVar(T, 0, 1))
        w_P1.append(M.addMVar(T, 0, 1))
        p_P1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P1[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump1 - sum(qpump1) == 0)
    M.addConstr(smp_P1    - sum(p_P1) == 0)
    
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    M.addConstr(IDD@Storage1 + A1@SumTurb1 - A1@SumPump1 - A1@SumTurb + A1@SumPump + A1@Spill1 == bStorage)
    UU = np.zeros(len(p_T1))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T1)):
        M.addConstr(IDD@u_T1[jj] - A1@v_T1[jj] + A1@w_T1[jj] == UU[jj])
        M.addConstr(IDD@u_P1[jj] - A1@v_P1[jj] + A1@w_P1[jj] == UU[jj])
        M.addConstr(IDD@u_T1[jj]+IDD@u_P1[jj] <= 1)
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P1[jj]-aj*IDD@qpump1[jj]          >= 0 )
        M.addConstr(IDD@p_T1[jj]-aj*IDD@q1[jj] - cj*IDD @ z1 <= bj)
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
        
        
        M.addConstr(IDD @ p_P1[jj] - Pmin01*IDD @ u_P1[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_P1[jj] - IDD @ p_P1[jj] >= 0)
        
        M.addConstr(IDD @ p_T1[jj] - Pmin01*IDD @ u_T1[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_T1[jj] - IDD @ p_T1[jj] >= 0)
    
        M.addConstr(IDD@u_T1[jj] - BU@v_T1[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T1[jj] + BT@w_T1[jj] <=onesVec-wAdd)
    
    M.addConstr(IDD @ Storage1-Capacity/2 * IDD @ z1          >= 0 )
    M.addConstr(IDD @ Storage1-Capacity/2 * IDD @ z1          <= Capacity/2 )
        
    
    #RHSc = np.sin(np.cumsum(np.ones(T))/20)*4000
    M.addConstr(IDD @ smp_T + IDD @ smp_T1 >= RHSc)
    M.addConstr(IDD @ smp_P + IDD @ smp_P1 >= 0.1*negativ(RHSc))
 
    
    
    
    positID = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    negatID = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    
    
    M.addConstr(IDD @ positID - positiv(RHSc) == 0)
    M.addConstr(IDD @ negatID - 0.1*negativ(RHSc) == 0)
    
    
    
        
    positDA = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    negatDA = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    
    M.addConstr(IDD @ positDA - IDD @ smp_T - IDD @ smp_T1 + positiv(RHSc) == 0)
    M.addConstr(IDD @ negatDA - IDD @ smp_P - IDD @ smp_P1 + 0.1*negativ(RHSc) == 0)
    
    
        
    OBJ = 0
    OBJ = OBJ + Frc @ positDA - Frc @ negatDA
    for sc in range(len(p_T)):
        OBJ = OBJ - Penlt*np.ones(T) @ v_T[sc] - Penlt*np.ones(T) @ w_T[sc]\
            - Penlt*np.ones(T)@v_P[sc] - Penlt*np.ones(T)@w_P[sc] 
    
    
    for sc in range(len(p_T1)):
        OBJ = OBJ - Penlt*np.ones(T)@v_T1[sc] - Penlt*np.ones(T)@w_T1[sc]\
            - Penlt*np.ones(T)@v_P1[sc] - Penlt*np.ones(T)@w_P1[sc] 
    OBJ = OBJ - Penlt*np.ones(T) @ Spill
    OBJ = OBJ - Penlt*np.ones(T) @ Spill1
    OBJ = OBJ + Frc2 @ positID - Frc2 @ negatID

    
    M.setParam('TimeLimit',TL)
    M.setObjective(OBJ, GRB.MAXIMIZE)
    
    
    
    for jj in range(len(p_T)):
        p_T[jj].start = ppT[jj]
        u_T[jj].start = uuT[jj]
        v_T[jj].start = vvT[jj]
        w_T[jj].start = wwT[jj]
        q[jj].start   = qq[jj]
    for jj in range(len(p_T1)):
        p_T1[jj].start = ppT1[jj]
        u_T1[jj].start = uuT1[jj]
        v_T1[jj].start = vvT1[jj]
        w_T1[jj].start = wwT1[jj]
        q1[jj].start   = qq1[jj]
        
    for jj in range(len(p_P)):
        p_P[jj].start = ppP[jj]
        u_P[jj].start = uuP[jj]
        v_P[jj].start = vvP[jj]
        w_P[jj].start = wwP[jj]
        qpump[jj].start   = qpm[jj]
    for jj in range(len(p_P1)):
        p_P1[jj].start = ppP1[jj]
        u_P1[jj].start = uuP1[jj]
        v_P1[jj].start = vvP1[jj]
        w_P1[jj].start = wwP1[jj]
        qpump1[jj].start   = qpm1[jj]
    
    Storage.start = Strg
    Spill.start   = Spll
    z.start       = zzzz
    SumTurb.start = SmTrb 
    SumPump.start = SmPmp
    
    Storage1.start = Strg1
    Spill1.start   = Spll1
    z1.start       = zzzz1
    SumTurb1.start = SmTrb1 
    SumPump1.start = SmPmp1
    
    M.optimize()

    return M.ObjVal, M



def OptiEntire(RHSc, assets, Frc1, Frc2, TL, env):
    Frc = Frc1
    T = len(Frc)
    L0 = 0
    PowerPl = 'HydroPowerPlant1'
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    Penlt = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['StartCost']
    
    Inflows = np.random.randint(2, size=len(Frc))
    Inflows = Inflows*Capacity/(180*sum(Inflows))
    bStorage = np.cumsum(Inflows*0)
    bStorage = L0 + np.array(bStorage)
    
    M = Model(env=env)
    Storage  = M.addMVar(T, 0, Capacity)
    Spill  = M.addMVar(T, 0, Capacity)
    z      = M.addMVar(T, 0, 1)
    z.vType = gp.GRB.BINARY
    
    SumPmax  = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    SmPowerMax1 = SumPmax1
    SumTurb = M.addMVar(T, 0, SumPmax)
    SumPump = M.addMVar(T, 0, SumPmax)
    
    smp_T   = M.addMVar(T, 0, SumPmax1)
    smp_P   = M.addMVar(T, 0, SumPmax1)
    
    
    
    
    q        = []
    p_T      = []
    u_T      = []
    v_T      = []
    w_T      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        q.append(M.addMVar(T, 0, assets['HydroPowerPlant1']['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T.append(M.addMVar(T, 0, 1))
        v_T.append(M.addMVar(T, 0, 1))
        w_T.append(M.addMVar(T, 0, 1))
        p_T.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumTurb - sum(q)   == 0)
    M.addConstr(smp_T   - sum(p_T) == 0)
    
    
    qpump = []
    p_P     = []
    u_P      = []
    v_P      = []
    w_P      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump.append(M.addMVar(T, 0, assets['HydroPowerPlant1']['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P.append(M.addMVar(T, 0, 1))
        v_P.append(M.addMVar(T, 0, 1))
        w_P.append(M.addMVar(T, 0, 1))
        p_P.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump - sum(qpump) == 0)
    M.addConstr(smp_P   - sum(p_P)   == 0)
    
    UT       = assets['HydroPowerPlant1']['TechnicalFeature']['Turbines']['Turbine1']['TminOn:']
    DT       = assets['HydroPowerPlant1']['TechnicalFeature']['Turbines']['Turbine1']['TminOff:']
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    M.addConstr(IDD@Storage + A1@SumTurb - A1@SumPump + A1@Spill == bStorage)
    UU = np.zeros(len(p_T))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T)):
        M.addConstr(IDD@u_T[jj] - A1@v_T[jj] + A1@w_T[jj] == UU[jj])
        M.addConstr(IDD@u_P[jj] - A1@v_P[jj] + A1@w_P[jj] == UU[jj])
        M.addConstr(IDD@u_T[jj]+IDD@u_P[jj] <= 1)
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P[jj]-aj*IDD@qpump[jj]          >= 0 )
        M.addConstr(IDD@p_T[jj]-aj*IDD@q[jj] - cj*IDD @ z <= bj)
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
        
        M.addConstr(IDD @ p_P[jj] - Pmin01*IDD @ u_P[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_P[jj] - IDD @ p_P[jj] >= 0)
        
        M.addConstr(IDD @ p_T[jj] - Pmin01*IDD @ u_T[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_T[jj] - IDD @ p_T[jj] >= 0)
    
        M.addConstr(IDD@u_T[jj] - BU@v_T[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T[jj] + BT@w_T[jj] <=onesVec-wAdd)
    
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          >= 0 )
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          <= Capacity/2 )
        
    
    
    PowerPl = 'HydroPowerPlant2'
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    
    Inflows = np.random.randint(2, size=len(Frc))
    Inflows = Inflows*Capacity/(180*sum(Inflows))
    bStorage = np.cumsum(Inflows*0)
    bStorage = L0 + np.array(bStorage)
    
    
    Storage1  = M.addMVar(T, 0, Capacity)
    Spill1  = M.addMVar(T, 0, Capacity)
    z1      = M.addMVar(T, 0, 1)
    z1.vType = gp.GRB.BINARY
    
    SumPmax = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    SmPowerMax2 = SumPmax1
    SumTurb1 = M.addMVar(T, 0, SumPmax)
    SumPump1 = M.addMVar(T, 0, SumPmax)
    smp_T1   = M.addMVar(T, 0, SumPmax1)
    smp_P1   = M.addMVar(T, 0, SumPmax1)
    
    q1        = []
    p_T1      = []
    u_T1      = []
    v_T1      = []
    w_T1      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        q1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T1.append(M.addMVar(T, 0, 1))
        v_T1.append(M.addMVar(T, 0, 1))
        w_T1.append(M.addMVar(T, 0, 1))
        p_T1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T1[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumTurb1 - sum(q1)   == 0)
    M.addConstr(smp_T1   - sum(p_T1) == 0)
    
    qpump1 = []
    p_P1     = []
    u_P1      = []
    v_P1      = []
    w_P1      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P1.append(M.addMVar(T, 0, 1))
        v_P1.append(M.addMVar(T, 0, 1))
        w_P1.append(M.addMVar(T, 0, 1))
        p_P1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P1[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump1 - sum(qpump1) == 0)
    M.addConstr(smp_P1    - sum(p_P1) == 0)
    
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    M.addConstr(IDD@Storage1 + A1@SumTurb1 - A1@SumPump1 - A1@SumTurb + A1@SumPump + A1@Spill1 == bStorage)
    UU = np.zeros(len(p_T1))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T1)):
        M.addConstr(IDD@u_T1[jj] - A1@v_T1[jj] + A1@w_T1[jj] == UU[jj])
        M.addConstr(IDD@u_P1[jj] - A1@v_P1[jj] + A1@w_P1[jj] == UU[jj])
        M.addConstr(IDD@u_T1[jj]+IDD@u_P1[jj] <= 1)
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P1[jj]-aj*IDD@qpump1[jj]          >= 0 )
        M.addConstr(IDD@p_T1[jj]-aj*IDD@q1[jj] - cj*IDD @ z1 <= bj)
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
        
        
        M.addConstr(IDD @ p_P1[jj] - Pmin01*IDD @ u_P1[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_P1[jj] - IDD @ p_P1[jj] >= 0)
        
        M.addConstr(IDD @ p_T1[jj] - Pmin01*IDD @ u_T1[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_T1[jj] - IDD @ p_T1[jj] >= 0)
    
        M.addConstr(IDD@u_T1[jj] - BU@v_T1[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T1[jj] + BT@w_T1[jj] <=onesVec-wAdd)
    
    M.addConstr(IDD @ Storage1-Capacity/2 * IDD @ z1          >= 0 )
    M.addConstr(IDD @ Storage1-Capacity/2 * IDD @ z1          <= Capacity/2 )
        
    
    #RHSc = np.sin(np.cumsum(np.ones(T))/20)*4000
    M.addConstr(IDD @ smp_T + IDD @ smp_T1 >= RHSc)
    M.addConstr(IDD @ smp_T + IDD @ smp_T1 >= 0.1*negativ(RHSc))

        
    positID = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    negatID = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    
    
    M.addConstr(IDD @ positID - positiv(RHSc) == 0)
    M.addConstr(IDD @ negatID - 0.1*negativ(RHSc) == 0)
    
    
    
        
    positDA = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    negatDA = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    
    M.addConstr(IDD @ positDA - IDD @ smp_T - IDD @ smp_T1 + positiv(RHSc) == 0)
    M.addConstr(IDD @ negatDA - IDD @ smp_P - IDD @ smp_P1 + 0.1*negativ(RHSc) == 0)
    
    
        
    OBJ = 0
    OBJ = OBJ + Frc @ positDA - Frc @ negatDA
    for sc in range(len(p_T)):
        OBJ = OBJ - Penlt*np.ones(T) @ v_T[sc] - Penlt*np.ones(T) @ w_T[sc]\
            - Penlt*np.ones(T)@v_P[sc] - Penlt*np.ones(T)@w_P[sc] 
    
    
    for sc in range(len(p_T1)):
        OBJ = OBJ - Penlt*np.ones(T)@v_T1[sc] - Penlt*np.ones(T)@w_T1[sc]\
            - Penlt*np.ones(T)@v_P1[sc] - Penlt*np.ones(T)@w_P1[sc] 
    OBJ = OBJ - Penlt*np.ones(T) @ Spill
    OBJ = OBJ - Penlt*np.ones(T) @ Spill1
    OBJ = OBJ + Frc2 @ positID - Frc2 @ negatID

    
    M.setParam('TimeLimit',TL)
    M.setObjective(OBJ, GRB.MAXIMIZE)
    
    
    
    M.optimize()

    return M.ObjVal, M

def transformAsset(assets, array):
    Keys1 = list(assets.keys())
    assets2 = copy.deepcopy(assets)
    for key in Keys1:
        Keys2 = list(assets[key]['TechnicalFeature']['Turbines'].keys())
        remainingArray = np.setdiff1d(range(len(Keys2)), array)
        for el in range(len(Keys2)):
            if el in remainingArray:
                del assets2[key]['TechnicalFeature']['Turbines'][Keys2[el]]  
    return assets2    





def BuildWS(ppT,    uuT,    vvT,    wwT,     qq,  \
        ppP,    uuP,    vvP,    wwP,     qpm, \
        ppT1,   uuT1,   vvT1,   wwT1,    qq1,  \
        ppP1,   uuP1,   vvP1,   wwP1,    qpm1, N1, N2, array):

    n1 = len(ppT)
    n2 = len(ppT1)
    AppT = [[]]*N1
    AuuT = [[]]*N1
    AvvT = [[]]*N1
    AwwT = [[]]*N1
    Aqq = [[]]*N1
    AppP = [[]]*N1
    AuuP = [[]]*N1
    AvvP = [[]]*N1
    AwwP = [[]]*N1
    Aqpm = [[]]*N1
    AppT1 = [[]]*N2
    AuuT1 = [[]]*N2
    AvvT1 = [[]]*N2
    AwwT1 = [[]]*N2
    Aqq1 = [[]]*N2
    AppP1 = [[]]*N2
    AuuP1 = [[]]*N2
    AvvP1 = [[]]*N2
    AwwP1 = [[]]*N2
    Aqpm1 = [[]]*N2  



    
    ii = 0
    for jj in range(N1):
        if jj in array:
            AppT[jj] = ppT[ii]
            AuuT[jj] = uuT[ii]
            AvvT[jj] = vvT[ii]
            AwwT[jj] = wwT[ii]
            Aqq[jj]  = qq[ii]
            AppP[jj] = ppP[ii]
            AuuP[jj] = uuP[ii]
            AvvP[jj] = vvP[ii]
            AwwP[jj] = wwP[ii]
            Aqpm[jj] = qpm[ii]

            ii = min(n1-1, ii + 1)
        else:
            AppT[jj] = ppT[ii]*0
            AuuT[jj] = uuT[ii]*0
            AvvT[jj] = vvT[ii]*0
            AwwT[jj] = wwT[ii]*0
            Aqq[jj]  = qq[ii]*0
            AppP[jj] = ppP[ii]*0
            AuuP[jj] = uuP[ii]*0
            AvvP[jj] = vvP[ii]*0
            AwwP[jj] = wwP[ii]*0
            Aqpm[jj] = qpm[ii]*0     

        
            
    ii = 0
    for jj in range(N2):
        if jj in array:
            AppT1[jj] = ppT1[ii]
            AuuT1[jj] = uuT1[ii]
            AvvT1[jj] = vvT1[ii]
            AwwT1[jj] = wwT1[ii]
            Aqq1[jj]  = qq1[ii]
            AppP1[jj] = ppP1[ii]
            AuuP1[jj] = uuP1[ii]
            AvvP1[jj] = vvP1[ii]
            AwwP1[jj] = wwP1[ii]
            Aqpm1[jj] = qpm1[ii]

            ii = min(n2-1, ii + 1)
        else:
            AppT1[jj] = ppT1[ii]*0
            AuuT1[jj] = uuT1[ii]*0
            AvvT1[jj] = vvT1[ii]*0
            AwwT1[jj] = wwT1[ii]*0
            Aqq1[jj]  = qq1[ii]*0
            AppP1[jj] = ppP1[ii]*0
            AuuP1[jj] = uuP1[ii]*0
            AvvP1[jj] = vvP1[ii]*0
            AwwP1[jj] = wwP1[ii]*0
            Aqpm1[jj] = qpm1[ii]*0   


    return AppT,    AuuT,    AvvT,    AwwT,     Aqq,  \
           AppP,    AuuP,    AvvP,    AwwP,     Aqpm, \
           AppT1,   AuuT1,   AvvT1,   AwwT1,    Aqq1,  \
           AppP1,   AuuP1,   AvvP1,   AwwP1,    Aqpm1

def BuildWS1(ppT,    uuT,    vvT,    wwT,     qq,  \
        ppP,    uuP,    vvP,    wwP,     qpm, \
        N1, array):

    n1 = len(ppT)

    AppT = [[]]*N1
    AuuT = [[]]*N1
    AvvT = [[]]*N1
    AwwT = [[]]*N1
    Aqq  = [[]]*N1
    AppP = [[]]*N1
    AuuP = [[]]*N1
    AvvP = [[]]*N1
    AwwP = [[]]*N1
    Aqpm = [[]]*N1



    
    ii = 0
    for jj in range(N1):
        if jj in array:
            AppT[jj] = ppT[ii]
            AuuT[jj] = uuT[ii]
            AvvT[jj] = vvT[ii]
            AwwT[jj] = wwT[ii]
            Aqq[jj]  = qq[ii]
            AppP[jj] = ppP[ii]
            AuuP[jj] = uuP[ii]
            AvvP[jj] = vvP[ii]
            AwwP[jj] = wwP[ii]
            Aqpm[jj] = qpm[ii]

            ii = min(n1-1, ii + 1)
        else:
            AppT[jj] = ppT[ii]*0
            AuuT[jj] = uuT[ii]*0
            AvvT[jj] = vvT[ii]*0
            AwwT[jj] = wwT[ii]*0
            Aqq[jj]  = qq[ii]*0
            AppP[jj] = ppP[ii]*0
            AuuP[jj] = uuP[ii]*0
            AvvP[jj] = vvP[ii]*0
            AwwP[jj] = wwP[ii]*0
            Aqpm[jj] = qpm[ii]*0     


    return AppT,    AuuT,    AvvT,    AwwT,     Aqq,  \
           AppP,    AuuP,    AvvP,    AwwP,     Aqpm

def OptimizeSinglePP(RHSc, assets, PowerPl, Frc1, Frc2, TL, env, Inflows, rhs3):
    Frc = Frc1
    T = len(Frc)
    L0 = 0
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    Penlt = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['StartCost']
    

    bStorage = np.cumsum(Inflows)
    bStorage = L0 + np.array(bStorage)
    
    M = Model(env=env)
    Storage  = M.addMVar(T, 0, Capacity)
    Spill  = M.addMVar(T, 0, Capacity)
    z      = M.addMVar(T, 0, 1)
    z.vType = gp.GRB.BINARY
    
    SumPmax = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    
    SumTurb = M.addMVar(T, 0, SumPmax)
    SumPump = M.addMVar(T, 0, SumPmax)
    
    smp_T  = M.addMVar(T, 0, SumPmax1)
    smp_P  = M.addMVar(T, 0, SumPmax1)
    
    q        = []
    p_T      = []
    u_T      = []
    v_T      = []
    w_T      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        q.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T.append(M.addMVar(T, 0, 1))
        v_T.append(M.addMVar(T, 0, 1))
        w_T.append(M.addMVar(T, 0, 1))
        p_T.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumTurb - sum(q)   == 0)
    M.addConstr(smp_T   - sum(p_T) == 0)
    
    qpump    = []
    p_P      = []
    u_P      = []
    v_P      = []
    w_P      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P.append(M.addMVar(T, 0, 1))
        v_P.append(M.addMVar(T, 0, 1))
        w_P.append(M.addMVar(T, 0, 1))
        p_P.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump - sum(qpump) == 0)
    M.addConstr(smp_P   - sum(p_P)   == 0)
    
    UT       = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['TminOn:']
    DT       = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['TminOff:']
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    M.addConstr(IDD@Storage + A1@SumTurb - A1@SumPump + A1@Spill == bStorage - rhs3)
    UU = np.zeros(len(p_T))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T)):
        M.addConstr(IDD@u_T[jj] - A1@v_T[jj] + A1@w_T[jj] == UU[jj])
        M.addConstr(IDD@u_P[jj] - A1@v_P[jj] + A1@w_P[jj] == UU[jj])
        M.addConstr(IDD@u_T[jj]+IDD@u_P[jj] <= 1)
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P[jj]-aj*IDD@qpump[jj]          >= 0 )
        M.addConstr(IDD@p_T[jj]-aj*IDD@q[jj] - cj*IDD @ z <= bj)
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
               
        M.addConstr(IDD @ p_P[jj] - Pmin01*IDD @ u_P[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_P[jj] - IDD @ p_P[jj] >= 0)
        
        M.addConstr(IDD @ p_T[jj] - Pmin01*IDD @ u_T[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_T[jj] - IDD @ p_T[jj] >= 0)
    
        M.addConstr(IDD@u_T[jj] - BU@v_T[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T[jj] + BT@w_T[jj] <=onesVec-wAdd)
    
    #RHSc = np.sin(np.cumsum(np.ones(T))/20)*4000
    M.addConstr(IDD @ smp_T   >= RHSc)
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          >= 0 )
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          <= Capacity/2 )
    
    
    M.addConstr(IDD @ smp_T >= RHSc)
    M.addConstr(IDD @ smp_T >= 0.1*negativ(RHSc))

        
    positID = M.addMVar(T, 0, SumPmax1 )
    negatID = M.addMVar(T, 0, SumPmax1 )
    
    
    M.addConstr(IDD @ positID - positiv(RHSc) == 0)
    M.addConstr(IDD @ negatID - 0.1*negativ(RHSc) == 0)

    


    positDA = M.addMVar(T, 0, SumPmax1)
    negatDA = M.addMVar(T, 0, SumPmax1)
    
    M.addConstr(IDD @ positDA - IDD @ smp_T  + positiv(RHSc) == 0)
    M.addConstr(IDD @ negatDA - IDD @ smp_P  + 0.1*negativ(RHSc) == 0)
    
    
        
    OBJ = 0
    OBJ = OBJ + Frc @ positDA - Frc @ negatDA
    for sc in range(len(p_T)):
        OBJ = OBJ - Penlt*np.ones(T) @ v_T[sc] - Penlt*np.ones(T) @ w_T[sc]\
            - Penlt*np.ones(T)@v_P[sc] - Penlt*np.ones(T)@w_P[sc] 
    
    
    OBJ = OBJ - Penlt*np.ones(T) @ Spill
    OBJ = OBJ + Frc2 @ positID - Frc2 @ negatID
    

        
    
    
    M.setParam('TimeLimit',TL)
    M.setObjective(OBJ, GRB.MAXIMIZE)
    M.optimize()
    ppT = []
    uuT = []
    vvT = []
    wwT = []
    qq  = []
    for jj in range(len(p_T)):
        ppT.append(p_T[jj].x)
        uuT.append(u_T[jj].x)
        vvT.append(v_T[jj].x)
        wwT.append(w_T[jj].x)
        qq.append(q[jj].x)
    ppP = []
    uuP = []
    vvP = []
    wwP = []
    qpm  = []
    for jj in range(len(p_P)):
        ppP.append(p_P[jj].x)
        uuP.append(u_P[jj].x)
        vvP.append(v_P[jj].x)
        wwP.append(w_P[jj].x)
        qpm.append(qpump[jj].x)
    Strg = Storage.x
    Spll = Spill.x
    zzzz = z.x
    SmTrb = SumTurb.x
    SmPmp = SumPump.x
    
    rhs2 = np.dot(A1, SmTrb ) - np.dot(A1,SmPmp)
  
    return Strg,  Spll,  zzzz,  SmTrb,  SmPmp,\
           ppT,   uuT,   vvT,   wwT,    qq, \
           ppP,   uuP,   vvP,   wwP,    qpm, M, rhs2




def OptimizeSinglePPWS(RHSc, assets, PowerPl, Frc1, Frc2, TL, env, Inflows, rhs3, Strg,  Spll,  zzzz,  SmTrb,  SmPmp, ppT,  uuT,   vvT,   wwT,    qq, ppP,   uuP,   vvP,   wwP,    qpm):
    Frc = Frc1
    T = len(Frc)
    L0 = 0
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    Penlt = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['StartCost']
    

    bStorage = np.cumsum(Inflows)
    bStorage = L0 + np.array(bStorage)
    
    M = Model(env=env)
    Storage  = M.addMVar(T, 0, Capacity)
    Spill  = M.addMVar(T, 0, Capacity)
    z      = M.addMVar(T, 0, 1)
    z.vType = gp.GRB.BINARY
    
    SumPmax = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    
    SumTurb = M.addMVar(T, 0, SumPmax)
    SumPump = M.addMVar(T, 0, SumPmax)
    
    smp_T  = M.addMVar(T, 0, SumPmax1)
    smp_P  = M.addMVar(T, 0, SumPmax1)
    
    q        = []
    p_T      = []
    u_T      = []
    v_T      = []
    w_T      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        q.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T.append(M.addMVar(T, 0, 1))
        v_T.append(M.addMVar(T, 0, 1))
        w_T.append(M.addMVar(T, 0, 1))
        p_T.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumTurb - sum(q) == 0)
    M.addConstr(smp_T   - sum(p_T) == 0)
    
    qpump = []
    p_P     = []
    u_P      = []
    v_P      = []
    w_P      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P.append(M.addMVar(T, 0, 1))
        v_P.append(M.addMVar(T, 0, 1))
        w_P.append(M.addMVar(T, 0, 1))
        p_P.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump - sum(qpump) == 0)
    M.addConstr(smp_P   - sum(p_P) == 0)
    
    UT       = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['TminOn:']
    DT       = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['TminOff:']
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    M.addConstr(IDD@Storage + A1@SumTurb - A1@SumPump + A1@Spill == bStorage - rhs3)
    UU = np.zeros(len(p_T))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T)):
        M.addConstr(IDD@u_T[jj] - A1@v_T[jj] + A1@w_T[jj] == UU[jj])
        M.addConstr(IDD@u_P[jj] - A1@v_P[jj] + A1@w_P[jj] == UU[jj])
        M.addConstr(IDD@u_T[jj]+IDD@u_P[jj] <= 1)
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P[jj]-aj*IDD@qpump[jj]          >= 0 )
        M.addConstr(IDD@p_T[jj]-aj*IDD@q[jj] - cj*IDD @ z <= bj)
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
               
        M.addConstr(IDD @ p_P[jj] - Pmin01*IDD @ u_P[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_P[jj] - IDD @ p_P[jj] >= 0)
        
        M.addConstr(IDD @ p_T[jj] - Pmin01*IDD @ u_T[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_T[jj] - IDD @ p_T[jj] >= 0)
    
        M.addConstr(IDD@u_T[jj] - BU@v_T[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T[jj] + BT@w_T[jj] <=onesVec-wAdd)
    
    #RHSc = np.sin(np.cumsum(np.ones(T))/20)*4000
    M.addConstr(IDD @ smp_T   >= RHSc)
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          >= 0 )
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          <= Capacity/2 )
    
    
    M.addConstr(IDD @ smp_T >= RHSc)
    M.addConstr(IDD @ smp_T >= 0.1*negativ(RHSc))

        
    positID = M.addMVar(T, 0, SumPmax1 )
    negatID = M.addMVar(T, 0, SumPmax1 )
    
    
    M.addConstr(IDD @ positID - positiv(RHSc) == 0)
    M.addConstr(IDD @ negatID - 0.1*negativ(RHSc) == 0)

    


    positDA = M.addMVar(T, 0, SumPmax1)
    negatDA = M.addMVar(T, 0, SumPmax1)
    
    M.addConstr(IDD @ positDA - IDD @ smp_T  + positiv(RHSc) == 0)
    M.addConstr(IDD @ negatDA - IDD @ smp_P  + 0.1*negativ(RHSc) == 0)
    
    
        
    OBJ = 0
    OBJ = OBJ + Frc @ positDA - Frc @ negatDA
    for sc in range(len(p_T)):
        OBJ = OBJ - Penlt*np.ones(T) @ v_T[sc] - Penlt*np.ones(T) @ w_T[sc]\
            - Penlt*np.ones(T)@v_P[sc] - Penlt*np.ones(T)@w_P[sc] 
    
    
    OBJ = OBJ - Penlt*np.ones(T) @ Spill
    OBJ = OBJ + Frc2 @ positID - Frc2 @ negatID
    

        
    
    
    M.setParam('TimeLimit',TL)
    M.setObjective(OBJ, GRB.MAXIMIZE)

    for jj in range(len(p_T)):
        p_T[jj].start = ppT[jj]
        u_T[jj].start = uuT[jj]
        v_T[jj].start = vvT[jj]
        w_T[jj].start = wwT[jj]
        q[jj].start   = qq[jj]
 
    for jj in range(len(p_P)):
        p_P[jj].start = ppP[jj]
        u_P[jj].start = uuP[jj]
        v_P[jj].start = vvP[jj]
        w_P[jj].start = wwP[jj]
        qpump[jj].start   = qpm[jj]
 
    
    Storage.start = Strg
    Spill.start   = Spll
    z.start       = zzzz
    SumTurb.start = SmTrb 
    SumPump.start = SmPmp
    
 
    
    
    
    M.optimize()
    ppT = []
    uuT = []
    vvT = []
    wwT = []
    qq  = []
    for jj in range(len(p_T)):
        ppT.append(p_T[jj].x)
        uuT.append(u_T[jj].x)
        vvT.append(v_T[jj].x)
        wwT.append(w_T[jj].x)
        qq.append(q[jj].x)
    ppP = []
    uuP = []
    vvP = []
    wwP = []
    qpm  = []
    for jj in range(len(p_P)):
        ppP.append(p_P[jj].x)
        uuP.append(u_P[jj].x)
        vvP.append(v_P[jj].x)
        wwP.append(w_P[jj].x)
        qpm.append(qpump[jj].x)
    Strg = Storage.x
    Spll = Spill.x
    zzzz = z.x
    SmTrb = SumTurb.x
    SmPmp = SumPump.x
    
    rhs2 = np.dot(A1, SmTrb ) - np.dot(A1,SmPmp)
  
    return Strg,  Spll,  zzzz,  SmTrb,  SmPmp,\
           ppT,   uuT,   vvT,   wwT,    qq, \
           ppP,   uuP,   vvP,   wwP,    qpm, M, rhs2







           
           
           




def OptiEntire3(RHSc, assets, Frc1, Frc2, TL, env, zzzz, u_T0, u_P0):
    Frc = Frc1
    T = len(Frc)
    L0 = 0
    PowerPl = 'HydroPowerPlant1'
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    Penlt = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['StartCost']
    
    Inflows = np.random.randint(2, size=len(Frc))
    Inflows = Inflows*Capacity/(180*sum(Inflows))
    bStorage = np.cumsum(Inflows*0)
    bStorage = L0 + np.array(bStorage)
    
    M = Model(env=env)
    Storage  = M.addMVar(T, 0, Capacity)
    Spill  = M.addMVar(T, 0, Capacity)
    z      = M.addMVar(T, 0, 1)
    #z.vType = gp.GRB.BINARY
    
    SumPmax  = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    SmPowerMax1 = SumPmax1
    SumTurb = M.addMVar(T, 0, SumPmax)
    SumPump = M.addMVar(T, 0, SumPmax)
    
    smp_T   = M.addMVar(T, 0, SumPmax1)
    smp_P   = M.addMVar(T, 0, SumPmax1)
    
    q        = []
    p_T      = []
    u_T      = []
    v_T      = []
    w_T      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        q.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T.append(M.addMVar(T, 0, 1))
        v_T.append(M.addMVar(T, 0, 1))
        w_T.append(M.addMVar(T, 0, 1))
        p_T.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        #u_T[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumTurb - sum(q) == 0)
    M.addConstr(smp_T   - sum(p_T) == 0)
    
    qpump = []
    p_P     = []
    u_P      = []
    v_P      = []
    w_P      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P.append(M.addMVar(T, 0, 1))
        v_P.append(M.addMVar(T, 0, 1))
        w_P.append(M.addMVar(T, 0, 1))
        p_P.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        #u_P[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump - sum(qpump) == 0)
    M.addConstr(smp_P   - sum(p_P) == 0)
    
    UT       = assets['HydroPowerPlant1']['TechnicalFeature']['Turbines']['Turbine1']['TminOn:']
    DT       = assets['HydroPowerPlant1']['TechnicalFeature']['Turbines']['Turbine1']['TminOff:']
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    M.addConstr(IDD@Storage + A1@SumTurb - A1@SumPump + A1@Spill == bStorage)
    UU = np.zeros(len(p_T))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T)):
        M.addConstr(IDD@u_T[jj] - A1@v_T[jj] + A1@w_T[jj] == UU[jj])
        M.addConstr(IDD@u_P[jj] - A1@v_P[jj] + A1@w_P[jj] == UU[jj])
        M.addConstr(IDD@u_T[jj] == u_T0[jj])
        M.addConstr(IDD@u_P[jj] == u_P0[jj])
        M.addConstr(IDD@u_T[jj]+IDD@u_P[jj] <= 1)
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P[jj]-aj*IDD@qpump[jj]          >= 0 )
        M.addConstr(IDD@p_T[jj]-aj*IDD@q[jj] - cj*IDD @ z <= bj)
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
        

        
        M.addConstr(IDD @ p_P[jj] - Pmin01*IDD @ u_P[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_P[jj] - IDD @ p_P[jj] >= 0)
        
        M.addConstr(IDD @ p_T[jj] - Pmin01*IDD @ u_T[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_T[jj] - IDD @ p_T[jj] >= 0)
    
        M.addConstr(IDD@u_T[jj] - BU@v_T[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T[jj] + BT@w_T[jj] <=onesVec-wAdd)
    
    M.addConstr(IDD@z == zzzz)
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          >= 0 )
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          <= Capacity/2 )    
     
    PowerPl = 'HydroPowerPlant2'
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    
    Inflows = np.random.randint(2, size=len(Frc))
    Inflows = Inflows*Capacity/(180*sum(Inflows))
    bStorage = np.cumsum(Inflows*0)
    bStorage = L0 + np.array(bStorage)
    
    
    Storage1  = M.addMVar(T, 0, Capacity)
    Spill1  = M.addMVar(T, 0, Capacity)
    z1      = M.addMVar(T, 0, 1)
    z1.vType = gp.GRB.BINARY
    
    SumPmax = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    SmPowerMax2 = SumPmax1
    SumTurb1 = M.addMVar(T, 0, SumPmax)
    SumPump1 = M.addMVar(T, 0, SumPmax)
    
    smp_T1 = M.addMVar(T, 0, SumPmax1)
    smp_P1 = M.addMVar(T, 0, SumPmax1)
    
    q1        = []
    p_T1      = []
    u_T1      = []
    v_T1      = []
    w_T1      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        q1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T1.append(M.addMVar(T, 0, 1))
        v_T1.append(M.addMVar(T, 0, 1))
        w_T1.append(M.addMVar(T, 0, 1))
        p_T1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T1[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumTurb1 - sum(q1)   == 0)
    M.addConstr(smp_T1   - sum(p_T1) == 0)
    
    qpump1 = []
    p_P1     = []
    u_P1      = []
    v_P1      = []
    w_P1      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P1.append(M.addMVar(T, 0, 1))
        v_P1.append(M.addMVar(T, 0, 1))
        w_P1.append(M.addMVar(T, 0, 1))
        p_P1.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P1[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump1 - sum(qpump1) == 0)
    M.addConstr(smp_P1   - sum(p_P1) == 0)
    
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    M.addConstr(IDD@Storage1 + A1@SumTurb1 - A1@SumPump1 - A1@SumTurb + A1@SumPump + A1@Spill1 == bStorage)
    UU = np.zeros(len(p_T1))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T1)):
        M.addConstr(IDD@u_T1[jj] - A1@v_T1[jj] + A1@w_T1[jj] == UU[jj])
        M.addConstr(IDD@u_P1[jj] - A1@v_P1[jj] + A1@w_P1[jj] == UU[jj])
        M.addConstr(IDD@u_T1[jj]+IDD@u_P1[jj] <= 1)
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P1[jj]-aj*IDD@qpump1[jj]          >= 0 )
        M.addConstr(IDD@p_T1[jj]-aj*IDD@q1[jj] - cj*IDD @ z1 <= bj)
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
        
        
        M.addConstr(IDD @ p_P1[jj] - Pmin01*IDD @ u_P1[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_P1[jj] - IDD @ p_P1[jj] >= 0)
        
        M.addConstr(IDD @ p_T1[jj] - Pmin01*IDD @ u_T1[jj] >= 0)
        M.addConstr(Pmax01*IDD @ u_T1[jj] - IDD @ p_T1[jj] >= 0)
    
        M.addConstr(IDD@u_T1[jj] - BU@v_T1[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T1[jj] + BT@w_T1[jj] <=onesVec-wAdd)
    
    #RHSc = np.sin(np.cumsum(np.ones(T))/20)*4000
    M.addConstr(IDD @ Storage1-Capacity/2 * IDD @ z1          >= 0 )
    M.addConstr(IDD @ Storage1-Capacity/2 * IDD @ z1          <= Capacity/2 )
    M.addConstr(IDD @ smp_T + IDD @ smp_T1 >= RHSc)
    M.addConstr(IDD @ smp_T + IDD @ smp_T1 >= 0.1*negativ(RHSc))

        
    positID = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    negatID = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    
    
    M.addConstr(IDD @ positID - positiv(RHSc) == 0)
    M.addConstr(IDD @ negatID - 0.1*negativ(RHSc) == 0)
    
    
    
        
    positDA = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    negatDA = M.addMVar(T, 0, SmPowerMax1+SmPowerMax2)
    
    M.addConstr(IDD @ positDA - IDD @ smp_T - IDD @ smp_T1 + positiv(RHSc) == 0)
    M.addConstr(IDD @ negatDA - IDD @ smp_P - IDD @ smp_P1 + 0.1*negativ(RHSc) == 0)
    
    
        
    OBJ = 0
    OBJ = OBJ + Frc @ positDA - Frc @ negatDA
    for sc in range(len(p_T)):
        OBJ = OBJ - Penlt*np.ones(T) @ v_T[sc] - Penlt*np.ones(T) @ w_T[sc]\
            - Penlt*np.ones(T)@v_P[sc] - Penlt*np.ones(T)@w_P[sc] 
    
    
    for sc in range(len(p_T1)):
        OBJ = OBJ - Penlt*np.ones(T)@v_T1[sc] - Penlt*np.ones(T)@w_T1[sc]\
            - Penlt*np.ones(T)@v_P1[sc] - Penlt*np.ones(T)@w_P1[sc] 
    OBJ = OBJ - Penlt*np.ones(T) @ Spill
    OBJ = OBJ - Penlt*np.ones(T) @ Spill1
    OBJ = OBJ + Frc2 @ positID - Frc2 @ negatID

    
    M.setParam('TimeLimit',TL)
    M.setObjective(OBJ, GRB.MAXIMIZE)
    
    
    
    M.optimize()

    ppT = []
    uuT = []
    vvT = []
    wwT = []
    qq  = []
    for jj in range(len(p_T)):
        ppT.append(p_T[jj].x)
        uuT.append(u_T[jj].x)
        vvT.append(v_T[jj].x)
        wwT.append(w_T[jj].x)
        qq.append(q[jj].x)
    ppP = []
    uuP = []
    vvP = []
    wwP = []
    qpm  = []
    for jj in range(len(p_P)):
        ppP.append(p_P[jj].x)
        uuP.append(u_P[jj].x)
        vvP.append(v_P[jj].x)
        wwP.append(w_P[jj].x)
        qpm.append(qpump[jj].x)
    Strg = Storage.x
    Spll = Spill.x
    zzzz = z.x
    SmTrb = SumTurb.x
    SmPmp = SumPump.x
    
    ppT1 = []
    uuT1 = []
    vvT1 = []
    wwT1 = []
    qq1  = []
    for jj in range(len(p_T1)):
        ppT1.append(p_T1[jj].x)
        uuT1.append(u_T1[jj].x)
        vvT1.append(v_T1[jj].x)
        wwT1.append(w_T1[jj].x)
        qq1.append(q1[jj].x)
    ppP1 = []
    uuP1 = []
    vvP1 = []
    wwP1 = []
    qpm1  = []
    for jj in range(len(p_P1)):
        ppP1.append(p_P1[jj].x)
        uuP1.append(u_P1[jj].x)
        vvP1.append(v_P1[jj].x)
        wwP1.append(w_P1[jj].x)
        qpm1.append(qpump1[jj].x)
    Strg1 = Storage1.x
    Spll1 = Spill1.x
    zzzz1 = z1.x
    SmTrb1 = SumTurb1.x
    SmPmp1 = SumPump1.x
    return Strg,  Spll,  zzzz,  SmTrb,  SmPmp,\
           ppT,   uuT,   vvT,   wwT,    qq, \
           ppP,   uuP,   vvP,   wwP,    qpm,\
           Strg1, Spll1, zzzz1, SmTrb1, SmPmp1,\
           ppT1,  uuT1,  vvT1,  wwT1,   qq1, \
           ppP1,   uuP1,  vvP1,  wwP1,   qpm1, M

      
def Visualize(SUM1, SPT, n1, n2, Y_labelSum, Y_labelPrice, TITLE, m):
    fig, ax1 = plt.subplots()
    n1 = 0
    n2 = len(SPT)
    
    color = 'tab:red'
    ax1.set_xlabel('Time (h)')
    ax1.set_ylabel(Y_labelSum, color=color)
    ax1.plot(SUM1[(n1-m):(n2-m)], color=color)
    ax1.tick_params(axis='y', labelcolor=color)
    ax1.set_ylim([0, 1.2*max(SUM1)])
    ax1.set_title(TITLE)
    
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    
    color = 'tab:blue'
    ax2.set_ylabel(Y_labelPrice, color=color)  # we already handled the x-label with ax1
    ax2.plot(SPT[n1:n2], color=color)
    ax2.tick_params(axis='y', labelcolor=color)
    
    fig.tight_layout()  # otherwise the right y-label is slightly clipped
    plt.show()     

def STR(val):
    if (val>100):
        return "%0.3f" % val
    else:
        return "%0.4f" % val
    
    
#The function OptiSingle writes out all the variables and equations associated with the hydro asset assets[PowerPl]
"""Description of the inputs
      RHSc    ... the Right-hand side of the coupling constraint: A*x >= RHSc
      assets  ... the list of assets. assets[PowerPl] is that single asset being optimized
      PowerPl ... the index the power plant being optimized
      Frc1    ... the hisorical time series of prices. Objective value coefficients
              ... Frc1 -> "1" is there because initially I planned to use ID prices Frc2. But this model contains only DA prices
      TL      ... time limit of optimuzation: M.setParam('TimeLimit',TL) 
      env     ... envvironment connection parameters
      Inflows ... the inflow vector (in the settings the inflow is zero: pump and release)
      RunOpti ... if RunOpti == 1, this function will add the coupling constraint and run the optimization: M.optimize()
      M       ... the inputted model M. M is gurobi model.. M.addMvar(), M.addConstr(),..., M.optimize()
      OBJ     ... the inputted objective function. It contains the prices and penalties for switches of the turbines, spills, and constraint valuations
"""
def OptiSingle(RHSc, assets, PowerPl, Frc1, TL, env, Inflows, RunOpti, M, OBJ):
    Frc = Frc1
    T = len(Frc)
    L0 = 0
    # the capacity of the power plant in cubic meters
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    Penlt = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['StartCost']
    
    # the RHS of the water balance equation
    bStorage = np.cumsum(Inflows)
    bStorage = L0 + np.array(bStorage)
    
    # M can be either some gurobi model M, or an empty variable []
    # if M is empty then it will be defined as follows
    if not M:
        M = Model(env=env)
    # otherwise M will be just updated
    
    # definition of the variables
    # variable STORAGE
    Storage  = M.addMVar(T, 0, Capacity)
    # variable SPILLAGE
    Spill  = M.addMVar(T, 0, Capacity)
    # Binary variable z: if STORAGE >= half-full, z=1, and z=0, otherwise
    z      = M.addMVar(T, 0, 1)
    z.vType = gp.GRB.BINARY
    
    
    SumPmax = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    
    # cumulative water release: m**3
    SumTurb = M.addMVar(T, 0, SumPmax)
    # cumulative water pump: m**3
    SumPump = M.addMVar(T, 0, SumPmax)
    
    # cumulave power due to the release in kWh
    smp_T  = M.addMVar(T, 0, SumPmax1)
    # cumulave power due to the pump in kWh
    smp_P  = M.addMVar(T, 0, SumPmax1)
    
    q        = []    # len(q) = the number of the turbines in the optimized asset: q[j][t] is the amount of water released through turbine j at time t
    p_T      = []    # the subscript T means the actions of turbines: p_T is the power produced by Turbine
    u_T      = []                                                   # u_T[j][t]=1 if turbine j is on at time t
    v_T      = []                                                   # u_T, v_T, and w_T refer to 3bin formulation 
    w_T      = []                                                   # u_T + u_P <= 1
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        q.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T.append(M.addMVar(T, 0, 1))
        v_T.append(M.addMVar(T, 0, 1))
        w_T.append(M.addMVar(T, 0, 1))
        p_T.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    # the sum of water releases
    M.addConstr(SumTurb - sum(q) == 0)
    # the sum of powers produced by all the turbines
    M.addConstr(smp_T   - sum(p_T) == 0)
    
    qpump = []
    p_P     = []
    u_P      = []
    v_P      = []
    w_P      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P.append(M.addMVar(T, 0, 1))
        v_P.append(M.addMVar(T, 0, 1))
        w_P.append(M.addMVar(T, 0, 1))
        p_P.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump - sum(qpump) == 0)
    M.addConstr(smp_P   - sum(p_P) == 0)
    
    
    
    UT       = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['TminOn']
    DT       = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['TminOff']
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    # Water balance equation
    M.addConstr(IDD@Storage + A1@SumTurb - A1@SumPump + A1@Spill == bStorage)
    UU = np.zeros(len(p_T))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T)):
        M.addConstr(IDD@u_T[jj] - A1@v_T[jj] + A1@w_T[jj] == UU[jj])    # u(t)-u(t-1) = v(t) - w(t)   for turbines
        M.addConstr(IDD@u_P[jj] - A1@v_P[jj] + A1@w_P[jj] == UU[jj])    # u(t)-u(t-1) = v(t) - w(t)   for pumps
        M.addConstr(IDD@u_T[jj]+IDD@u_P[jj] <= 1)                       # the turbine can either release or pump, but not both
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P[jj]-aj*IDD@qpump[jj]          >= 0 )        #Pump:    power = a(j)*Q....  
        M.addConstr(IDD@p_T[jj]-aj*IDD@q[jj] - cj*IDD @ z <= bj)        #Release: power = a(j)*Q + c(j)*z + b(j). z is binary var: if storage level >= Capacity/2, then z = 1, else z = 0
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
               
        M.addConstr(IDD @ p_P[jj] - Pmin01*IDD @ u_P[jj] >= 0)    # p - Pmin*u >=0
        M.addConstr(Pmax01*IDD @ u_P[jj] - IDD @ p_P[jj] >= 0)    # Pmax*u - p >=0
        
        M.addConstr(IDD @ p_T[jj] - Pmin01*IDD @ u_T[jj] >= 0)    # p - Pmin*u >=0
        M.addConstr(Pmax01*IDD @ u_T[jj] - IDD @ p_T[jj] >= 0)    # Pmax*u - p >=0
        # Timing Constraints
        # UT       = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['TminOn']
        # DT       = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['TminOff']
        # A2       = np.tril(np.ones([T,T]),-UT)
        # BU       = A1 - A2 
        # A2       = np.tril(np.ones([T,T]),-DT)
        # BT       = A1 - A2    
        M.addConstr(IDD@u_T[jj] - BU@v_T[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T[jj] + BT@w_T[jj] <=onesVec-wAdd)
    
    # the following two constraints are equivalent to the condition:
    # if STORAGE >= half-full, z=1, and z=0, otherwise
   
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          >= 0 )
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          <= Capacity/2 )
    # if RunOpti = 1, then the coupling constraint is inputted and the optimization is launched:
    if RunOpti:
        M.addConstr(IDD @ smp_T >= positiv(RHSc))
        M.addConstr(IDD @ smp_P >= negativ(RHSc))

    # introduction of the objective function
    if not OBJ: 
        OBJ = 0
    # adding objective value coefficients:
    OBJ = OBJ + Frc @ smp_T - Frc @ smp_P
    # adding penalties for switches
    for sc in range(len(p_T)):
        OBJ = OBJ - Penlt*np.ones(T) @ v_T[sc] - Penlt*np.ones(T) @ w_T[sc]\
            - Penlt*np.ones(T)@v_P[sc] - Penlt*np.ones(T)@w_P[sc] 
    
    # adding penalties for spill
    OBJ = OBJ - Penlt*np.ones(T) @ Spill
    
    # stting the time limit
    M.setParam('TimeLimit',TL)
    # final setting up the model
    M.setObjective(OBJ, GRB.MAXIMIZE)
    # Running the model
    if RunOpti:
        M.optimize()

    # this function returns: the model, the objective, the storage, z, total sum power (turbines), total sum power (pumps)
    # total sum water (turbines), total sum water (pumps), Spillage.
    return M, OBJ, Storage, z, smp_T, smp_P, SumTurb, SumPump, Spill





















#The function OptiSingle writes out all the variables and equations associated with the hydro asset assets[PowerPl]
"""Description of the inputs
      RHSc    ... the Right-hand side of the coupling constraint: A*x >= RHSc
      assets  ... the list of assets. assets[PowerPl] is that single asset being optimized
      PowerPl ... the index the power plant being optimized
      Frc1    ... the hisorical time series of prices. Objective value coefficients
              ... Frc1 -> "1" is there because initially I planned to use ID prices Frc2. But this model contains only DA prices
      TL      ... time limit of optimuzation: M.setParam('TimeLimit',TL) 
      env     ... envvironment connection parameters
      Inflows ... the inflow vector (in the settings the inflow is zero: pump and release)
      RunOpti ... if RunOpti == 1, this function will add the coupling constraint and run the optimization: M.optimize()
      M       ... the inputted model M. M is gurobi model.. M.addMvar(), M.addConstr(),..., M.optimize()
      OBJ     ... the inputted objective function. It contains the prices and penalties for switches of the turbines, spills, and constraint valuations
"""
def OptiSingle1(RHSc, assets, PowerPl, Frc1, TL, env, Inflows, RunOpti, M, OBJ):
    Frc = Frc1
    T = len(Frc)
    L0 = 0
    # the capacity of the power plant in cubic meters
    Capacity = assets[PowerPl]['TechnicalFeature']['Capacity']
    Penlt = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['StartCost']
    
    # the RHS of the water balance equation
    bStorage = np.cumsum(Inflows)
    bStorage = L0 + np.array(bStorage)
    
    # M can be either some gurobi model M, or an empty variable []
    # if M is empty then it will be defined as follows
    if not M:
        M = Model(env=env)
    # otherwise M will be just updated
    
    # definition of the variables
    # variable STORAGE
    Storage  = M.addMVar(T, 0, Capacity)
    # variable SPILLAGE
    Spill  = M.addMVar(T, 0, Capacity)
    # Binary variable z: if STORAGE >= half-full, z=1, and z=0, otherwise
    z      = M.addMVar(T, 0, 1)
    z.vType = gp.GRB.BINARY
    
    
    SumPmax = 0
    SumPmax1 = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines'].keys():
        SumPmax  = SumPmax  + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']
        SumPmax1 = SumPmax1 + assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']
    
    # cumulative water release: m**3
    SumTurb = M.addMVar(T, 0, SumPmax)
    # cumulative water pump: m**3
    SumPump = M.addMVar(T, 0, SumPmax)
    
    # cumulave power due to the release in kWh
    smp_T  = M.addMVar(T, 0, SumPmax1)
    # cumulave power due to the pump in kWh
    smp_P  = M.addMVar(T, 0, SumPmax1)
    
    q        = []    # len(q) = the number of the turbines in the optimized asset: q[j][t] is the amount of water released through turbine j at time t
    p_T      = []    # the subscript T means the actions of turbines: p_T is the power produced by Turbine
    u_T      = []                                                   # u_T[j][t]=1 if turbine j is on at time t
    v_T      = []                                                   # u_T, v_T, and w_T refer to 3bin formulation 
    w_T      = []                                                   # u_T + u_P <= 1
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        q.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_T.append(M.addMVar(T, 0, 1))
        v_T.append(M.addMVar(T, 0, 1))
        w_T.append(M.addMVar(T, 0, 1))
        p_T.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_T[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    # the sum of water releases
    M.addConstr(SumTurb - sum(q) == 0)
    # the sum of powers produced by all the turbines
    M.addConstr(smp_T   - sum(p_T) == 0)
    
    qpump = []
    p_P     = []
    u_P      = []
    v_P      = []
    w_P      = []
    jj = 0
    for tb in assets[PowerPl]['TechnicalFeature']['Turbines']:
        qpump.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Xbar']))
        u_P.append(M.addMVar(T, 0, 1))
        v_P.append(M.addMVar(T, 0, 1))
        w_P.append(M.addMVar(T, 0, 1))
        p_P.append(M.addMVar(T, 0, assets[PowerPl]['TechnicalFeature']['Turbines'][tb]['Pmax']))
        u_P[jj].vType = gp.GRB.BINARY
        jj = jj + 1
    M.addConstr(SumPump - sum(qpump) == 0)
    M.addConstr(smp_P   - sum(p_P) == 0)
    
    
    
    UT       = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['TminOn']
    DT       = assets[PowerPl]['TechnicalFeature']['Turbines']['Turbine1']['TminOff']
    IDD      = np.eye(T)
    A1       = np.tril(np.ones([T,T]),0)
    A2       = np.tril(np.ones([T,T]),-UT)
    BU       = A1 - A2 
    A2       = np.tril(np.ones([T,T]),-DT)
    BT       = A1 - A2  
    
    zeroVec         = np.array([0]*T)
    onesVec         = np.array([1]*T)
    vAdd = np.array(CUMSUMV(np.zeros(DT), T))
    wAdd = np.array(CUMSUMV(np.zeros(DT), T))
    
    # Water balance equation
    M.addConstr(IDD@Storage + A1@SumTurb - A1@SumPump + A1@Spill == bStorage)
    UU = np.zeros(len(p_T))
    KEY = list(assets[PowerPl]['TechnicalFeature']['Turbines'].keys())
    for jj in range(len(p_T)):
        M.addConstr(IDD@u_T[jj] - A1@v_T[jj] + A1@w_T[jj] == UU[jj])    # u(t)-u(t-1) = v(t) - w(t)   for turbines
        M.addConstr(IDD@u_P[jj] - A1@v_P[jj] + A1@w_P[jj] == UU[jj])    # u(t)-u(t-1) = v(t) - w(t)   for pumps
        M.addConstr(IDD@u_T[jj]+IDD@u_P[jj] <= 1)                       # the turbine can either release or pump, but not both
        aj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['ajT']
        bj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['bjT']
        cj = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['cjT']
        M.addConstr(IDD@p_P[jj]-aj*IDD@qpump[jj]          >= 0 )        #Pump:    power = a(j)*Q....  
        M.addConstr(IDD@p_T[jj]-aj*IDD@q[jj] - cj*IDD @ z <= bj)        #Release: power = a(j)*Q + c(j)*z + b(j). z is binary var: if storage level >= Capacity/2, then z = 1, else z = 0
        Pmin01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmin']
        Pmax01 = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['Pmax']
               
        M.addConstr(IDD @ p_P[jj] - Pmin01*IDD @ u_P[jj] >= 0)    # p - Pmin*u >=0
        M.addConstr(Pmax01*IDD @ u_P[jj] - IDD @ p_P[jj] >= 0)    # Pmax*u - p >=0
        
        M.addConstr(IDD @ p_T[jj] - Pmin01*IDD @ u_T[jj] >= 0)    # p - Pmin*u >=0
        M.addConstr(Pmax01*IDD @ u_T[jj] - IDD @ p_T[jj] >= 0)    # Pmax*u - p >=0
        # Timing Constraints
        UT       = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['TminOn']
        DT       = assets[PowerPl]['TechnicalFeature']['Turbines'][KEY[jj]]['TminOff']
        A2       = np.tril(np.ones([T,T]),-UT)
        BU       = A1 - A2 
        A2       = np.tril(np.ones([T,T]),-DT)
        BT       = A1 - A2    
        M.addConstr(IDD@u_T[jj] - BU@v_T[jj] >=zeroVec+vAdd)
        M.addConstr(IDD@u_T[jj] + BT@w_T[jj] <=onesVec-wAdd)
    
    # the following two constraints are equivalent to the condition:
    # if STORAGE >= half-full, z=1, and z=0, otherwise
   
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          >= 0 )
    M.addConstr(IDD @ Storage-Capacity/2 * IDD @ z          <= Capacity/2 )
    # if RunOpti = 1, then the coupling constraint is inputted and the optimization is launched:
    if RunOpti:
        M.addConstr(IDD @ smp_T >= positiv(RHSc))
        M.addConstr(IDD @ smp_P >= negativ(RHSc))

    # introduction of the objective function
    if not OBJ: 
        OBJ = 0
    # adding objective value coefficients:
    OBJ = OBJ + Frc @ smp_T - Frc @ smp_P
    # adding penalties for switches
    for sc in range(len(p_T)):
        OBJ = OBJ - Penlt*np.ones(T) @ v_T[sc] - Penlt*np.ones(T) @ w_T[sc]\
            - Penlt*np.ones(T)@v_P[sc] - Penlt*np.ones(T)@w_P[sc] 
    
    # adding penalties for spill
    OBJ = OBJ - Penlt*np.ones(T) @ Spill
    
    # stting the time limit
    M.setParam('OutputFlag', 0)
    M.setParam('TimeLimit',TL)
    # final setting up the model
    M.setObjective(OBJ, GRB.MAXIMIZE)
    # Running the model
    if RunOpti:
        M.optimize()

    # this function returns: the model, the objective, the storage, z, total sum power (turbines), total sum power (pumps)
    # total sum water (turbines), total sum water (pumps), Spillage.
    return M, OBJ, Storage, z, smp_T, smp_P, SumTurb, SumPump, Spill, p_T, p_P, u_T, u_P, v_T, v_P, w_T, w_P




def transform(p_T):
    m = len(p_T)
    OutPut = []
    for j in range(m):
        OutPut.append(p_T[j].x)
    return OutPut


def Define_MILP_model(assets, RHSc, Frc, kn, env):
    KEYS = list(assets.keys())
    STORAGE = []; Z = []; SMP_T = []; SMP_P = []; SUMTURB = []; SUMPUMP = []; SPILLAGE = []; UT = []; UP = []
    #M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0                                                 = OptiSingle( RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0, p_T0, p_P0, u_T0, u_P0, v_T0, v_P0, w_T0, w_P0  = OptiSingle1(RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])
    T = len(Frc)
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

def Define_MILP_model0(assets, RHSc, Frc, kn, env):
    KEYS = list(assets.keys())
    STORAGE = []; Z = []; SMP_T = []; SMP_P = []; SUMTURB = []; SUMPUMP = []; SPILLAGE = []; UT = []; UP = []
    #M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0                                                 = OptiSingle( RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])
    M, OBJ, Storage0, z0, smp_T0, smp_P0, SumTurb0, SumPump0, Spill0, p_T0, p_P0, u_T0, u_P0, v_T0, v_P0, w_T0, w_P0  = OptiSingle1(RHSc, assets, KEYS[0], Frc, 60, env, RHSc*0, 0, [], [])
    T = len(Frc)
    #XP  = M.addMVar(T, 0, 10**8)
    #XT  = M.addMVar(T, 0, 10**8)
    #OBJ = OBJ - np.ones(T)*10000@XP  - np.ones(T)*10000@XT
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
    M.addConstr(sum(SMP_T) >= positiv(RHSc))
    M.addConstr(sum(SMP_P) >= negativ(RHSc))
    return M, STORAGE, Z, SMP_T, SMP_P, SUMTURB, SUMPUMP, SPILLAGE, UT, UP

def Gradual_Increase0(assets, RHSc, Frc, TimeLimits, SubPools, u_TA, u_PA, env):
    M=[[]]*len(TimeLimits); STORAGE=[[]]*len(TimeLimits); Z=[[]]*len(TimeLimits); SMP_T=[[]]*len(TimeLimits); SMP_P=[[]]*len(TimeLimits); SUMTURB=[[]]*len(TimeLimits); SUMPUMP=[[]]*len(TimeLimits); SPILLAGE=[[]]*len(TimeLimits); UT=[[]]*len(TimeLimits); UP=[[]]*len(TimeLimits)
    m = len(TimeLimits)
    ii = 0
    M0, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model0(assets, RHSc, Frc, SubPools[ii],env)
    M[ii] = M0; Z[ii] = Z0; STORAGE[ii] = STORAGE0; SMP_T[ii] = SMP_T0; SMP_P[ii] = SMP_P0
    SUMTURB[ii] = SUMTURB0; SUMPUMP[ii] = SUMPUMP0; SPILLAGE[ii] = SPILLAGE0; UT[ii] = UT0; UP[ii] = UP0
    M[ii].setParam('TimeLimit',TimeLimits[ii])
    M[ii].setParam('OutputFlag', 1)
    M[ii].optimize() 
    for ii in range(1, m):
        M0, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model0(assets, RHSc, Frc, SubPools[ii],env)
        M[ii] = M0; Z[ii] = Z0; STORAGE[ii] = STORAGE0; SMP_T[ii] = SMP_T0; SMP_P[ii] = SMP_P0
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
        M[ii].setParam('OutputFlag', 1)
        M[ii].optimize() 
    return M[-1]

# def Gradual_Increase(assets, RHSc, Frc, TimeLimits, SubPools, u_TA, u_PA, env):
#     M=[[]]*len(TimeLimits); STORAGE=[[]]*len(TimeLimits); Z=[[]]*len(TimeLimits); SMP_T=[[]]*len(TimeLimits); SMP_P=[[]]*len(TimeLimits); SUMTURB=[[]]*len(TimeLimits); SUMPUMP=[[]]*len(TimeLimits); SPILLAGE=[[]]*len(TimeLimits); UT=[[]]*len(TimeLimits); UP=[[]]*len(TimeLimits)
#     m = len(TimeLimits)
#     ii = 0
#     M0, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model(assets, RHSc, Frc, SubPools[ii],env)
#     M[ii] = M0; Z[ii] = Z0; STORAGE[ii] = STORAGE0; SMP_T[ii] = SMP_T0; SMP_P[ii] = SMP_P0
#     SUMTURB[ii] = SUMTURB0; SUMPUMP[ii] = SUMPUMP0; SPILLAGE[ii] = SPILLAGE0; UT[ii] = UT0; UP[ii] = UP0
#     M[ii].setParam('TimeLimit',TimeLimits[ii])
#     M[ii].setParam('OutputFlag', 1)
#     M[ii].optimize() 
#     for ii in range(1, m):
#         M0, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0, UT0, UP0 = Define_MILP_model(assets, RHSc, Frc, SubPools[ii],env)
#         M[ii] = M0; Z[ii] = Z0; STORAGE[ii] = STORAGE0; SMP_T[ii] = SMP_T0; SMP_P[ii] = SMP_P0
#         SUMTURB[ii] = SUMTURB0; SUMPUMP[ii] = SUMPUMP0; SPILLAGE[ii] = SPILLAGE0; UT[ii] = UT0; UP[ii] = UP0
        
#         for j in range(SubPools[ii-1]):
#             Z[ii][j].start = Z[ii-1][j].x
#             m = len(UT[ii-1][j])
#             for i in range(m):
#                 UT[ii][j][i].start = UT[ii-1][j][i].x
#                 UP[ii][j][i].start = UP[ii-1][j][i].x
                
        
#         for j in range(SubPools[ii-1], SubPools[ii]):
#             m = len(UT[ii][j])
#             for i in range(m):
#                 UT[ii][j][i].start = u_TA[j][i]
#                 UP[ii][j][i].start = u_PA[j][i]    
#         M[ii].setParam('TimeLimit',TimeLimits[ii])
#         M[ii].setParam('OutputFlag', 1)
#         M[ii].optimize() 
#     return M[-1], SMP_T[-1], SMP_P[-1], Z[-1], STORAGE[-1], UT[-1], UP[-1]

def Gradual_Increase(assets, RHSc, Frc, TimeLimits, SubPools, u_TA, u_PA, env):
    M=[[]]*len(TimeLimits); STORAGE=[[]]*len(TimeLimits); Z=[[]]*len(TimeLimits); 
    SMP_T=[[]]*len(TimeLimits); SMP_P=[[]]*len(TimeLimits); SUMTURB=[[]]*len(TimeLimits); 
    SUMPUMP=[[]]*len(TimeLimits); SPILLAGE=[[]]*len(TimeLimits); 
    UT=[[]]*len(TimeLimits); UP=[[]]*len(TimeLimits)
    m = len(TimeLimits)
    ii = 0
    M0, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0,\
        UT0, UP0 = Define_MILP_model(assets, RHSc, Frc, SubPools[ii],env)
    M[ii] = M0; Z[ii] = Z0; STORAGE[ii] = STORAGE0; SMP_T[ii] = SMP_T0; 
    SMP_P[ii] = SMP_P0;     SUMTURB[ii] = SUMTURB0; SUMPUMP[ii] = SUMPUMP0; 
    SPILLAGE[ii] = SPILLAGE0; UT[ii] = UT0; UP[ii] = UP0
    M[ii].setParam('TimeLimit',TimeLimits[ii])
    M[ii].setParam('OutputFlag', 1)
    M[ii].optimize() 
    for ii in range(1, m):
        M0, STORAGE0, Z0, SMP_T0, SMP_P0, SUMTURB0, SUMPUMP0, SPILLAGE0,\
            UT0, UP0 = Define_MILP_model(assets, RHSc, Frc, SubPools[ii],env)
        M[ii] = M0; Z[ii] = Z0; STORAGE[ii] = STORAGE0; 
        SMP_T[ii] = SMP_T0; SMP_P[ii] = SMP_P0
        SUMTURB[ii] = SUMTURB0; SUMPUMP[ii] = SUMPUMP0; 
        SPILLAGE[ii] = SPILLAGE0; UT[ii] = UT0; UP[ii] = UP0
        
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
        M[ii].setParam('OutputFlag', 1)
        M[ii].optimize() 
    return M[-1], SMP_T[-1], SMP_P[-1], Z[-1], STORAGE[-1], UT[-1], UP[-1]

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
