# -*- coding: utf-8 -*-
"""
Created on Thu Feb 29 14:45:08 2024

@author: brgris
"""
import numpy as np
import argparse
#import argparse

def indices_linha(lin, nome):
    origem = lin[:,0]     # Nome da Barra de origem
    destino = lin[:,1]   # Nome da Barra de destino
    ifr = np.zeros(len(lin), dtype=np.int16) # Indice da barra de origem
    ito = np.zeros(len(lin), dtype=np.int16) # Indice da barra de destino
    # Determinação das barras de origem e de destino
    for cont in range(len(lin)):
        for k in range(len(nome)):
            if nome[k]== origem[cont]:
                ifr[cont]=k
            if nome[k]== destino[cont]:
                ito[cont]=k
    return (ifr, ito)

def matriz_inc(indice_origem, indice_destino, n_bar):
    n_linha = len(indice_origem)    
    Ainc = np.zeros((n_linha,n_bar), dtype=np.int16) #Matriz de incidência
    Af = np.zeros((n_linha,n_bar), dtype=np.int16) #Matriz de incidência das barras de origem
    At = np.zeros((n_linha,n_bar), dtype=np.int16) #Matriz de incidência das barras de destino
    for k in range(n_linha):
        Ainc[k,indice_origem[k]]=1
        Ainc[k,indice_destino[k]]=-1
        Af[k,indice_origem[k]]= 1
        At[k,indice_destino[k]]= 1
    return Ainc, Af, At

def matriz_admitancia (lin, barra, Ainc):
    shcar = np.array(barra[:,10], dtype=np.float64) # Carga reativa na barra [pu]
    R = np.array(lin[:,3], dtype=np.float64) #Resistência*100%
    R = R/100 #Resistência [pu]
    X = np.array(lin[:,4], dtype=np.float64) #Reatância*100%
    X = X/100 #Reatância [pu]
    Ysh_2 = np.array(lin[:,5], dtype=np.float64) #y_shunt total*100%
    Ysh_2 = Ysh_2/200 #y_shunt em um dos terminais [pu]
    Ybus = np.zeros((len(barra),len(barra)), dtype=complex)
    Yp = np.zeros((len(lin),len(lin)), dtype=complex)
    Yshunt = np.zeros((len(lin)), dtype=complex)
    Yp = np.linalg.inv(np.diag(R+1j*X)) #Obtenção da matriz de admitância primitiva
    Yshunt = 1j*Ysh_2  #Yshunt
    Ybus = (Ainc.T)@(Yp)@(Ainc) + np.diag(np.abs(Ainc.T)@(Yshunt)) + np.diag(1j*shcar)  # Obtenção da matriz de admitância Ybus    
    # print(Ybus)
    return Ybus, Yp

def parametros(V, ang, gkm, bkm, bar_index, barPQ, Pcalc, Qcalc):
    n = len(V)
    H = np.zeros((n,n), dtype=float)
    for k in range(n):
        if tipo[k] !=2:
            H[k,k]=-Qcalc[k]-bkm[k,k]*V[k]**2
            for m in range(n):
                delta = ang[k]-ang[m]
                if k != m:
                    H[k,m]=V[k]*V[m]*(gkm[k,m]*np.sin(delta)-bkm[k,m]*np.cos(delta))
     
    L = np.zeros((n,n), dtype=float)
    for k in range(n):
        if tipo[k] == 0:
            L[k,k]= Qcalc[k]/V[k]-bkm[k,k]*V[k]
            for m in range(n):
                delta = ang[k]-ang[m]
                if k != m:
                    L[k,m]=V[k]*(gkm[k,m]*np.sin(delta)-bkm[k,m]*np.cos(delta))
                    
    M = np.zeros((n,n), dtype=float)
    for k in range(n):
        if tipo[k] == 0:
            M[k,k]= Pcalc[k]-gkm[k,k]*V[k]**2
            for m in range(n):
                delta = ang[k]-ang[m]
                if k != m:
                    M[k,m]=-V[k]*V[m]*(gkm[k,m]*np.cos(delta)+bkm[k,m]*np.sin(delta))
                    
    N = np.zeros((n,n), dtype=float)
    for k in range(n):
        if tipo[k] !=2:
            N[k,k]= Pcalc[k]/V[k]+gkm[k,k]*V[k]
            for m in range(n):
                delta = ang[k]-ang[m]
                if k != m:
                    N[k,m]=V[k]*(gkm[k,m]*np.cos(delta)+bkm[k,m]*np.sin(delta)) 
    H = H[np.ix_(bar_index,bar_index)]
    L = L[np.ix_(barPQ,barPQ)]
    M = M[np.ix_(barPQ,bar_index)]
    N = N[np.ix_(bar_index,barPQ)]                
    return H, L, M, N

def calculo_residuo(V, gkm, bkm, ang, tipo, barra):
    pcar = np.array(barra[:,8], dtype=np.float64)
    pcar = pcar/base # Potência ativa demandada [pu]
    qcar = np.array(barra[:,9], dtype=np.float64)
    qcar = qcar/base # Potência reativa demandada [pu]
    n = len(barra)
    Pcalc = np.zeros(n)
    Qcalc = np.zeros(n)
    for k in range(n):
        if tipo[k] != 2:
           for m in range(n):
               delta = ang[k]-ang[m]
               Pcalc[k]=Pcalc[k]+V[k]*V[m]*(gkm[k,m]*np.cos(delta)+bkm[k,m]*np.sin(delta))
               Qcalc[k]=Qcalc[k]+V[k]*V[m]*(gkm[k,m]*(np.sin(delta))-(bkm[k,m]*np.cos(delta)))
               dP = pger-pcar-Pcalc #Cálculo do resíduo de potência ativa
               dQ = qger-qcar-Qcalc #Cálculo do resíduo de potência reativa
               dP[k]= pger[k]-pcar[k]-Pcalc[k] 
               dQ[k]= qger[k]-qcar[k]-Qcalc[k]
    # max_dP = max(abs(dP))   
    # max_dQ = max(abs(dQ)) 
    return dP, dQ, Pcalc, Qcalc

# def metodo_desacoplado(H,L,M,N,bar_index,barPQ):
def calculo_ang_desacoplado(H, bar_index, ang, dP):    
    dP = dP[np.ix_(bar_index)]
    # dQ[np.ix_(barPQ)] = dQ[np.ix_(barPQ)]
    ddelta = np.linalg.solve(H,dP.T)
    maxdP = max(abs(dP))  
    ang[np.ix_(bar_index)] = ang[np.ix_(bar_index)]+ddelta
    return ang, maxdP


def calculo_V_desacoplado(L, barPQ, V, dQ):    
    dQ = dQ[np.ix_(barPQ)]
    dV = np.linalg.solve(L,dQ.T)
    maxdQ = max(abs(dQ))  
    V[np.ix_(barPQ)]=V[np.ix_(barPQ)]+dV 
    return V, maxdQ


def define_indices(barra):
    NPV=0
    NPQ=0
    for cont in range(len(barra)):
        if tipo[cont]==1:
            NPV += 1
        if tipo[cont]==0:    
            NPQ += 1
    
    bar_index = np.zeros((NPQ+NPV)) # Índice de barras PQ ou PV
    bar_index = np.array(bar_index,dtype=np.int16)
    
    barPV = []
    barPQ = []
    bar_index=[]
    
    for cont in range(len(barra)):
        if tipo[cont]==1:
            barPV.append(cont)
            bar_index.append(cont)
        if tipo[cont]==0:
            barPQ.append(cont)
            bar_index.append(cont)
    return NPV, NPQ, bar_index, barPV, barPQ

def fluxo_potencia_newton (barra, Ybus, max_iter=20):
    gkm = Ybus.real
    bkm = Ybus.imag
    maxdPQ=1      
    num_iter=0
    V = np.array(barra[:,2], dtype=np.float64) # Modulo da V nas barras [pu] - Estimativa inicial
    ang = np.array(barra[:,3], dtype=np.float64) # Angulo inicial das barras
    NPV, NPQ, bar_index, barPV, barPQ = define_indices(barra)
    while (num_iter < max_iter) and (maxdPQ>0.003):
        dP, dQ, Pcalc, Qcalc = calculo_residuo(V, gkm, bkm, ang, tipo, barra)
        H, L, M, N = parametros(V, ang, gkm, bkm, bar_index, barPQ, Pcalc, Qcalc)
        J = np.block([[H, N],[M, L]])  
        dPQ = np.block([dP[np.ix_(bar_index)], dQ[np.ix_(barPQ)]])
        maxdPQ = max(abs(dPQ))      
        dX = np.linalg.solve(J,dPQ.T)   
        ddelta = dX[0:(NPQ+NPV)]
        ang[np.ix_(bar_index)] = ang[np.ix_(bar_index)]+ddelta
        if NPQ>0:
            dV = dX[(NPQ+NPV):(2*NPQ+NPV)]
            V[np.ix_(barPQ)]=V[np.ix_(barPQ)]+dV   
        print(f'Iter_P: {num_iter} - Erro: {maxdPQ:.5f}')
        num_iter+=1 
        for i in range(len(barra)):  
            print(f'Barra: {barra[i,1]} - Tensão: {V[i]:6.3f} - Angulo:{ang[i]:6.3f}')
            
            
def fluxo_potencia_desacoplado(barra, Ybus, max_iter=30):
    gkm = Ybus.real
    bkm = Ybus.imag
    maxdP=1 
    maxdQ=1
    num_iter_P=0
    num_iter_Q=0
    num_iter=0
    V = np.array(barra[:,2], dtype=np.float64) # Modulo da V nas barras [pu] - Estimativa inicial
    ang = np.array(barra[:,3], dtype=np.float64) # Angulo inicial das barras
    NPV, NPQ, bar_index, barPV, barPQ = define_indices(barra)
    
    while (num_iter < max_iter):
        # dP, dQ, Pcalc, Qcalc = calculo_residuo(V, gkm, bkm, ang, tipo, barra)
        if maxdP>0.003:
            dP, dQ, Pcalc, Qcalc = calculo_residuo(V, gkm, bkm, ang, tipo, barra)
            H, L, M, N = parametros(V, ang, gkm, bkm, bar_index, barPQ, Pcalc, Qcalc)            
            ang, maxdP = calculo_ang_desacoplado(H, bar_index, ang, dP)
            print(f'Iter_P: {num_iter_P} - Erro: {maxdP:.5f}') 
            num_iter_P += 1
            for i in range(len(barra)):  
                print(f'Barra: {barra[i,1]} - Tensão: {V[i]:6.3f} - Angulo:{ang[i]:6.3f}')
        if maxdQ>0.003:
            dP, dQ, Pcalc, Qcalc = calculo_residuo(V, gkm, bkm, ang, tipo, barra)
            H, L, M, N = parametros(V, ang, gkm, bkm, bar_index, barPQ, Pcalc, Qcalc)
            V, maxdQ = calculo_V_desacoplado(L, barPQ, V, dQ)
            print(f'Iter_Q: {num_iter_Q} - Erro: {maxdQ:.5f}')  
            num_iter_Q += 1
            for i in range(len(barra)):  
                print(f'Barra: {barra[i,1]} - Tensão: {V[i]:6.3f} - Angulo:{ang[i]:6.3f}')
        num_iter+=1
            
def fluxo_potencia_desacoplado_rapido(barra, Ybus, max_iter=30):  
    
    tipo = np.array(barra[:,0], dtype = np.int16) # Tipo das barras
    gkm = Ybus.real
    bkm = Ybus.imag
    maxdP=1 
    maxdQ=1
    num_iter_P=0
    num_iter_Q=0
    num_iter=0
    V = np.array(barra[:,2], dtype=np.float64) # Modulo da V nas barras [pu] - Estimativa inicial
    ang = np.array(barra[:,3], dtype=np.float64) # Angulo inicial das barras
    NPV, NPQ, bar_index, barPV, barPQ = define_indices(barra)
    
    while (num_iter < max_iter):
        # dP, dQ, Pcalc, Qcalc = calculo_residuo(V, gkm, bkm, ang, tipo, barra)
        if maxdP>0.003:
            dP, dQ, Pcalc, Qcalc = calculo_residuo(V, gkm, bkm, ang, tipo, barra)
            dP=dP/V
            ang, maxdP = calculo_ang_desacoplado(-bkm[np.ix_(bar_index,bar_index)], bar_index, ang, dP)
            print(f'Iter_P: {num_iter_P} - Erro: {maxdP:.5f}') 
            num_iter_P += 1
            for i in range(len(barra)):  
                print(f'Barra: {barra[i,1]} - Tensão: {V[i]:6.3f} - Angulo:{ang[i]:6.3f}')
        if maxdQ>0.003:
            dP, dQ, Pcalc, Qcalc = calculo_residuo(V, gkm, bkm, ang, tipo, barra)
            dQ=dQ/V
            V, maxdQ = calculo_V_desacoplado(-bkm[np.ix_(barPQ,barPQ)], barPQ, V, dQ)
            print(f'Iter_Q: {num_iter_Q} - Erro: {maxdQ:.5f}')  
            num_iter_Q += 1
            for i in range(len(barra)):  
                print(f'Barra: {barra[i,1]} - Tensão: {V[i]:6.3f} - Angulo:{ang[i]:6.3f}')
        num_iter+=1  
        
#                 Tipo Nome   V[pu] Ang[rad] Pg[MW] Qg[Mvar] Qmin[Mvar] Qmax[Mvar] Pl[MW] Ql[Mvar] Sh[Mvar]   Vbase
barra = np.array([[2, 'B1',   1.00,    0.,   0.000,   0.,   -1000.,      1000.,     0.000,   0.000,     0,       230],
                  [0, 'B2',   1.00,    0.,   0.000,   0.,   -1000.,      1000.,     800.0,   280.0,     0,       230],
                  [1, 'B3',   1.05,    0.,   520.0,   0.,   -300.,       336.,      80.00,   40.00,     0,       230],
                  [0, 'B4',   1.00,    0.,   0.000,   0.,   -1000.,      1000.,     0.000,   0.000,     0,       230],
                  [0, 'B5',   1.02,    0.,   0.000,   0.,   -1000.,      1000.,     0.000,   0.000,     0,       230]])


#                   De   Para           R%        X%       Shunt
linha = np.array([['B1', 'B5',  1. ,   0.150 ,   2.00   ,  0.000 ],
                  ['B3', 'B4',  1. ,   0.075 ,   1.00   ,  0.000 ],
                  ['B2', 'B4',  1. ,   0.900 ,   10.0   ,  172.0 ],
                  ['B2', 'B5',  1. ,   0.450 ,   5.0    ,  88.00 ],
                  ['B4', 'B5',  1. ,   0.225 ,   2.5    ,  44.00]])

#                 Tipo Nome   V[pu] Ang[rad] Pg[MW] Qg[Mvar] Qmin[Mvar] Qmax[Mvar] Pl[MW] Ql[Mvar] Sh[Mvar]   Vbase
# barra = np.array([[2, 'B1',   1.00,    0.,   0.000,   0.,   -1000.,      1000.,     0.00,     0.00,     0,       230],
#                   [0, 'B2',   1.00,    0.,   0.000,   0.,   -1000.,      1000.,     30.00,   -7.00,     0,       230]])

# #                  De   Para           R%        X%       Shunt
# linha = np.array([['B1', 'B2',  1. ,   20.0 ,   100.0  ,  4.00 ]])


def main(metodo):
    global base, tipo, pger, qger
    base=100
    lin = np.array(linha)
    ifr, ito = indices_linha(lin, barra[:,1])
    Ainc, Af, At = matriz_inc(ifr, ito, len(barra)) # Determinação das Matrizes de incidência
    Ybus, Yp = matriz_admitancia (lin, barra, Ainc) 

    tipo = np.array(barra[:,0], dtype = np.int16) # Tipo das barras
    pger = np.array(barra[:,4], dtype=np.float64) # Potência ativa gerada nas barras (definidas pelo usuário)
    pger = pger/base # Potência ativa gerada [pu]
    qger = np.array(barra[:,5], dtype=np.float64) # Potência reativa gerada nas barras (definidas pelo usuário)
    qger = qger/base # Potência reativa gerada [pu]
    qmin = np.array(barra[:,6], dtype=np.float64)
    qmin = qmin/base # Potência reativa mínima que pode ser gerada [pu]
    qmax = np.array(barra[:,7], dtype=np.float64)
    qmax = qmax/base # Potência reativa máxima que pode ser gerada [pu]
    pcar = np.array(barra[:,8], dtype=np.float64)
    pcar = pcar/base # Potência ativa demandada [pu]
    qcar = np.array(barra[:,9], dtype=np.float64)
    qcar = qcar/base # Potência reativa demandada [pu]
    
    if metodo=='NR':
       fluxo_potencia_newton (barra, Ybus, max_iter=20)   
    if metodo=='D':
       fluxo_potencia_desacoplado(barra, Ybus, max_iter=30)
    if metodo=='DR':
       fluxo_potencia_desacoplado_rapido(barra, Ybus, max_iter=30) 
        
# main('NR') 
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Escolha um método para execução.")
    parser.add_argument("metodo", choices=["NR", "D", "DR"], help="Método a ser utilizado (NR, D, ou DR)")
    args = parser.parse_args()
    main(args.metodo)    
    
    
# !python metodos_FP.py NR    
# !python metodos_FP.py D    
# !python metodos_FP.py DR        
    
    
    
    