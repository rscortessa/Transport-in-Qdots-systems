import numpy as np
import re
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm
# Se lee el archivo aqui
def read_file(filename,filename2,n):
    archivo=open(filename)
    #print(archivo.read())
    Data=[[] for i in range(n)]
    A=archivo.read().split("\n")
    #print(A)
    for line in A:
        C=re.findall("-*"+"[0-9]+.[0-9]+",line)
        if len(C)==n:
            for j in range(n):
                Data[j].append(float(C[j]))
    try:
        archiva=open(filename2)
        #print(archiva.read())
        A=archiva.read().split("\n")
        for line in A:
            C=re.findall("-*"+"[0-9]+.[0-9]+",line)
            if len(C)==n:
                for j in range(n):
                    Data[j].append(float(C[j]))
        return Data
    except:     
        return Data
# Realiza una regresión lineal.
def LinearR2(col_1, col_2):
    col_1 = sm.add_constant(col_1)
    model = sm.OLS(col_2, col_1)
    results = model.fit()
    COEF = [results.params[0], results.bse[0], results.params[1], results.bse[1],results.rsquared]
    COEF =  np.array(COEF)
    return COEF
def colores():
    return ["red","blue","orange","yellow","gray","green","magenta","black","cyan","purple","beige", "blueviolet","brick","darkorange","gold3","orchid4","tan","aqua"]
### AQUI SE CREAN LAS FUNCIONES PARA EL CÁLCULO DE LOS MÁXIMOS EN LA FIGURA DE MÉRITO
def Maximos_ZT(ZT,E,mu):
    ZT_Aux=[ZT[i] for i in range(int(len(E)/2))]
    E_aux=[E[i] for i in range(int(len(E)/2))]
    indice=ZT_Aux.index(max(ZT_Aux))
    return [E_aux[indice],mu-(E_aux[indice]-mu),ZT_Aux[indice]]

def all_Maximos_ZT(T,datos,mu):
    maximos=np.zeros((4,len(T)))
    for i in range(len(T)):
        maximos[0,i]=T[i]
        maximos[1,i]=Maximos_ZT(datos[i][8],datos[i][0],mu)[0]
        maximos[2,i]=Maximos_ZT(datos[i][8],datos[i][0],mu)[1]
        maximos[3,i]=Maximos_ZT(datos[i][8],datos[i][0],mu)[2]
    return maximos

def density_ZT(T,datos):
    m=len(T)
    n=len(datos[0][0])
    k=len(datos[11][0])
    print(m,n,k)
    density=np.zeros((3,m*n))
    for i in range(m):
        for j in range(n):
            density[0,i*n+j]=T[i]
            density[1,i*n+j]=datos[i][0][j]
            density[2,i*n+j]=datos[i][8][j]
    return density
