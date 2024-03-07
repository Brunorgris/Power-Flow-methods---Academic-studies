
"""

ALGORITMO PARA CÁLCULO DE FLUXO DE POTÊNCIA
AUTOR: BRUNO RAFAEL GRIS
   
"""
import numpy as np
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

def define_indices(barra):
    tipo = np.array(barra[:,0], dtype = np.int16) # Tipo das barras
    classe = np.array(barra[:,12], dtype = str) # Tipo das barras
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
    barGer = []

    for cont in range(len(barra)):
        if tipo[cont]==1:
            barPV.append(cont)
            bar_index.append(cont)
        if tipo[cont]==0:
            barPQ.append(cont)
            bar_index.append(cont)
            
    for cont in range(len(barra)):
        if classe[cont]=='G':
            barGer.append(cont)
  
                
    return NPV, NPQ, bar_index, barPV, barPQ, barGer

def calc_potencia(tensao, ang, Ybus):
        Vc = tensao*np.exp(1j*ang) #Cálculo do vetor de tensão complexa
        Scalc = Vc*(np.conj(Ybus@Vc)) #Cálculo da potência complexa injetada nas barras
        Pcalc = np.real(Scalc) #Cálculo da potência ativa injetada nas barras
        Qcalc = np.imag(Scalc) #Cálculo da potência reativa injetada nas barras
        return Vc, Scalc, Pcalc, Qcalc   

def fluxo_potencia (barra, Ybus, base=100, max_iter=20):
    global qmax, qmin, barPQ, barPV, qcar
    """
    INICIALIZAÇÃO DE VARIÁVEIS
    """
    tensao = np.array(barra[:,2], dtype=np.float64) # Modulo da tensao nas barras [pu] - Estimativa inicial
    ang = np.array(barra[:,3], dtype=np.float64) # Angulo inicial das barras
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
    
    NPV, NPQ, bar_index, barPV, barPQ, barGer = define_indices(barra)
            
    #Inicialização de variáveis do processo iterativo
    maxdPQ=1
    num_iter = 0
    # Inicio do processo iterativo
    while (num_iter < max_iter) and (maxdPQ>0.001):
        num_iter = num_iter+1  #incremento do contador de numero de iterações
        Vc, Scalc, Pcalc, Qcalc = calc_potencia(tensao, ang, Ybus)
        qger = np.array(barra[:,5], dtype=np.float64) # Potência reativa gerada nas barras (definidas pelo usuário)
        qger = qger/base # Potência reativa gerada [pu]
        qger_calc = Qcalc+qcar # Q injetado é igual Qcalc e portando, igual a qger-qcar. Logo o qger tem que ser qinj+qcar
        for i in range(len(barGer)):
             if qger_calc[barGer[i]] < qmin[barGer[i]] or qger_calc[barGer[i]] > qmax[barGer[i]]:
                 print('Houve violação do limite reativo na barra:', barra[barGer[i],1])
                 barra[barGer[i],0]=0
                 if qger_calc[barGer[i]] < qmin[barGer[i]]:
                     barra[barGer[i],5]=barra[barGer[i],6]
                 if qger_calc[barGer[i]] > qmax[barGer[i]]:
                     barra[barGer[i],5]=barra[barGer[i],7]
             if qger_calc[barGer[i]] > qmin[barGer[i]] and qger_calc[barGer[i]] < qmax[barGer[i]]:
                 barra[barGer[i],0]=1
                 tensao[barGer[i]]=barra[barGer[i],2]
             
        Vc, Scalc, Pcalc, Qcalc = calc_potencia(tensao, ang, Ybus)     
        NPV, NPQ, bar_index, barPV, barPQ, barGer = define_indices(barra)  
        dP = pger-pcar-Pcalc #Cálculo do resíduo de potência ativa
        dQ = qger-qcar-Qcalc #Cálculo do resíduo de potência reativa
        #Obtenção da matriz Jacobiana
        Dva = 1j*np.diag(Vc) #Cálculo da derivada da tensão complexa em função do ângulo
        Dvv = np.diag(np.exp(1j*ang)) #Cálculo da derivada da tensão complexa em função da tensão em módulo
        Dsa = np.diag(Vc)@(np.conj(Ybus@Dva))+np.diag(np.conj(Ybus@Vc))@Dva #Cálculo da derivada da Potência complexa em função do ângulo
        Dsv = np.diag(Vc)@(np.conj(Ybus@Dvv))+np.diag(np.conj(Ybus@Vc))@Dvv #Cálculo da derivada da Potência complexa em função da tensão
        Dpa = np.real(Dsa) #Obtenção da matriz de derivada da Potência ativa em função do ângulo
        Dpv = np.real(Dsv) #Obtenção da matriz de derivada da Potência ativa em função do tensão
        Dqa = np.imag(Dsa) #Obtenção da matriz de derivada da Potência reativa em função do ângulo
        Dqv = np.imag(Dsv) #Obtenção da matriz de derivada da Potência reativa em função do tensão
        Dpa = Dpa[np.ix_(bar_index,bar_index)] #Indexação da matriz de derivada parcial Dpa
        Dqv = Dqv[np.ix_(barPQ,barPQ)] #Indexação da matriz de derivada parcial Dpv
        Dqa = Dqa[np.ix_(barPQ,bar_index)] #Indexação da matriz de derivada parcial Dqa
        Dpv = Dpv[np.ix_(bar_index,barPQ)] #Indexação da matriz de derivada parcial Dqv
        J = np.block([[Dpa, Dpv],[Dqa, Dqv]]) #Matriz Jacobiana
        #Obtenção da direção para novo processo iterativo 
        dP[np.ix_(bar_index)] = dP[np.ix_(bar_index)] #Indexação do vetor de resíduos de potência ativa
        # print(', '.join(map(str, np.round(dP,5))))
        dQ[np.ix_(barPQ)] = dQ[np.ix_(barPQ)] #Indexação do vetor de resíduos de potência reativa
        dPQ = np.block([dP[np.ix_(bar_index)], dQ[np.ix_(barPQ)]]) 
        maxdPQ = max(abs(dPQ)) #Norma infinito do vetor de resíduos de potência ativa e reativa
        dX = np.linalg.solve(J,dPQ.T)
        #Atualização das variáveis
        ddelta = dX[0:(NPQ+NPV)]
        ang[np.ix_(bar_index)] = ang[np.ix_(bar_index)]+ddelta
        if NPQ>0:
            dV = dX[(NPQ+NPV):(2*NPQ+NPV)]
            tensao[np.ix_(barPQ)]=tensao[np.ix_(barPQ)]+dV
            tensao = abs(tensao)
    #Término do processo iterativo e obtenção das grandezas de interesse
    return tensao, ang    


#def fluxpot(disjuntores, args):
def main():
    base = 100 # Potencia base do sistema
    #Dados barra - 0-Barra PQ, 1- Barra PV, 2 - Barra de referência
    
    #                 Tipo Nome   V[pu] Ang[rad] Pg[MW] Qg[Mvar] Qmin[Mvar] Qmax[Mvar] Pl[MW] Ql[Mvar] Sh[Mvar]   Vbase
    # barra = np.array([[2, 'B1',   1.00,    0.,   0.000,   0.,   -1000.,      1000.,     0.000,   0.000,     0,       230],
    #                   [0, 'B2',   1.00,    0.,   0.000,   0.,   -1000.,      1000.,     800.0,   280.0,     0,       230],
    #                   [1, 'B3',   1.05,    0.,   520.0,   0.,   -300.,       336.,      80.00,   40.00,     0,       230],
    #                   [0, 'B4',   1.00,    0.,   0.000,   0.,   -1000.,      1000.,     0.000,   0.000,     0,       230],
    #                   [0, 'B5',   1.02,    0.,   0.000,   0.,   -1000.,      1000.,     0.000,   0.000,     0,       230]])



    # #                   De   Para           R%        X%       Shunt
    # linha = np.array([['B1', 'B5',  1. ,   0.150 ,   2.00   ,  0.000 ],
    #                   ['B3', 'B4',  1. ,   0.075 ,   1.00   ,  0.000 ],
    #                   ['B2', 'B4',  1. ,   0.900 ,   10.0   ,  172.0 ],
    #                   ['B2', 'B5',  1. ,   0.450 ,   5.0    ,  88.00 ],
    #                   ['B4', 'B5',  1. ,   0.225 ,   2.5    ,  44.00]])
    
    barra = np.array([[2, 'B1',   1.00,    0.,   0.000,   0.,   -9999.,      9999.,     0.000,   0.000,     0,       230,   'C'],
                      [0, 'B2',   1.00,    0.,   0.000,   0.,   -1000.,      1000.,     800.0,   280.0,     0,       230,   'C'],
                      [1, 'B3',   1.05,    0.,   520.0,   0.,   -300,        340.,      80.00,   40.00,     0,       230,   'G'],
                      [0, 'B4',   1.00,    0.,   0.000,   0.,   -1000.,      1000.,     0.000,   0.000,     0,       230,   'C'],
                      [0, 'B5',   1.02,    0.,   0.000,   0.,   -1000.,      1000.,     0.000,   0.000,     0,       230,   'C']])


    #                   De   Para           R%        X%       Shunt
    linha = np.array([['B1', 'B5',  1. ,   0.150 ,   2.00   ,  0.000 ],
                      ['B3', 'B4',  1. ,   0.075 ,   1.00   ,  0.000 ],
                      ['B2', 'B4',  1. ,   0.900 ,   10.0   ,  172.0 ],
                      ['B2', 'B5',  1. ,   0.450 ,   5.0    ,  88.00 ],
                      ['B4', 'B5',  1. ,   0.225 ,   2.5    ,  44.00]])
    
    lin = np.array(linha)
    ifr, ito = indices_linha(lin, barra[:,1])
    Ainc, Af, At = matriz_inc(ifr, ito, len(barra)) # Determinação das Matrizes de incidência
    Ybus, Yp = matriz_admitancia (lin, barra, Ainc) 
    tensao, ang = fluxo_potencia(barra, Ybus)
    # Vbase = np.array(barra[:,11], dtype=np.float64) #tensão base
    Vc = tensao*np.exp(1j*ang) #Cálculo do vetor de tensão complexa
    Scalc = Vc*(np.conj(Ybus@Vc)) #Cálculo da potência complexa injetada nas barras
    Pcalc = np.real(Scalc) #Cálculo da potência ativa injetada nas barras
    Qcalc = np.imag(Scalc) #Cálculo da potência reativa injetada nas barras
    # Ikm = abs(Yp@Ainc@Vc*base*1000/(np.sqrt(3)*Af@Vbase))
    # Imk = -Yp@Ainc@Vc*base*1000/(np.sqrt(3)*At@Vbase)
    Skm = (Af@Vc)*np.conj(Yp@Ainc@Vc)
    Pkm = np.real(Skm)
    Qkm = np.imag(Skm)
    Pinj = np.round(Pcalc*base,0)
    Qinj = np.round(Qcalc*base,0)
    # Ikm=np.round(Ikm,0)
    Pkm=np.round(Pkm*base,0)
    Qkm=np.round(Qkm*base,0)
    # tensao = np.round(tensao*Vbase,2)
    nomebar = barra[:,1]
    return tensao, ang, nomebar, barra, Pinj, Qinj, Pkm, Qkm

if __name__== "__main__":
    tensao, ang, nomebar, barra, Pinj, Qinj, Pkm, Qkm = main()
    for i in range(len(barra)):  
        print(f'Barra: {nomebar[i]} - Tensão: {tensao[i]:6.3f} - Angulo:{ang[i]:6.3f}')
