from numpy import *
from math import *
#Ne pas oublier l'initialisation , module azote dispo azote dispo
def fmin(n,T) :
    if T[n]==0 :
        return 0
    else :
        s=0
        m=len(T)
        for i in range m:
            if T[i]>0 :
                s+=T[i]
        return T[m]/s
def Tmoy(T) :
    m=len(T)
    s=0
    n=0
    for i in range m:
        if T[i]>0
        s+=T[i]
        n+=1
    return s/n
        
#fonction temps normalisé
def Mresidu(Mr,fmin):
    s=Mr*fmin
    return s
#residu azoté disponible
def Mcompoorga(Ma,fmin):
    s=Ma*fmin
    return s
#composé organique
def Mhumus(fdsc,a,Tmoy,Norga,Da,T) :
    if T<!0 :
        return 0
    else :
        Mh=fdsc*?*exp(a*(Tmoy-Tref)*Norga*Da*ep*1000*(Beta/((%arg+?)*(CaCO3+d)))