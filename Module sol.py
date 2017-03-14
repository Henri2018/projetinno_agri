from numpy import *
from math import *
#Ne pas oublier l'initialisation , module azote dispo azote dispo
def fmin(n,T) :
#T la temperature et n le jour
    if T[n]==0 :
        return 0
    else :
        s=0
        m=len(T)
        for i in range m:
            if T[i]>0 :
                s+=T[i]
        return T[m]/s
def Tmoy(T,n) :
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
#residu azoté disponible, Mr est le coefficient de mineralisation des residus de culture
def Mcompoorga(Ma,fmin):
    s=Ma*fmin
    return s
#composé organique, Ma coefficient de mineralisation des composes organiques
def Mhumus(fdsc,a,beta,delta,Tmoy,Tref,Norga,Da,ep,Beta,pcarg,CaC03,d) :
    #pc arg représente le % d'argile
    if T<!0 :
        return 0
    else :
        Mh=fdsc*beta*exp(a*(Tmoy-Tref)*Norga*Da*ep*1000*(Beta/((pcarg+delta)*(CaCO3+d))))
        return Mh
#fdc le coefficient relatif à la gestion de la parcelle considérée
#Norganique représente le pourcentage d’azote total de la parcelle
#Da la densité apparente de l’horizon labouré 
# ep l’épaisseur de l’horizon labouré 
# %arg et CaCO3, respectivement les pourcentages d’argile et de calcaire de l’horizon labouré

def CaU(mu,MS,MS7,SDT,SDT1,l) :
    #l et mu représentent respectivement l’ordonnée à l’origine et la pente de la relation linéaire entre le CAU et la vitesse de croissance de la culture dans les sept jours précédents l’apport [Ln(MSj/10) – ln(MSj-7/10)/(SDTj – SDTj-7)].
    CaU=min(100,l+(mu*100*(((ln(MS/10) – ln(MS7/10))*(MS/10))/(SDT–SDT7))))
    return CaU