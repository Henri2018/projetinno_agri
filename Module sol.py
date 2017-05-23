import numpy as np
import math as math
#Ne pas oublier l'initialisation , module azote dispo azote dispo
def fmin(n,T) : #modifier pour prendre une température limite as forcemment egale à 0, regarder les bibliotheques +faire tests avec valeurs speciales
#T la temperature et n le jour
    if T[n]==0 :
        return 0
    else :
        return np.average(T)/n
        
#fonction temps normalisé 
def Mresidu(fmin):
    s=0.016*fmin
    return s
#residu azoté disponible, Mr est le coefficient de mineralisation des residus de culture http://www.supagro.fr/ress-pepites/matiereorganique/co/3_EntreesMO.html
def Mcompoorga(fmin):
    s=0.02*fmin
    return s
#composé organique, Ma coefficient de mineralisation des composes organiques http://www.groupe-frayssinet.fr/Actualites/Mineralisation-de-la-matiere-organique
def Mhumus(T) :
    #pc arg représente le % d'argile
    if np.average(T)<0 :
        return 0
    else :
        Mh=0.23*math.exp(0.115*(np.average(T)-15))*1000*(1.3/(11*60))
        return Mh
    
def sol(n,T,X):
    f=fmin(n,T)
    T1=np.average(T)
    Sol=Mresidu(f)+Mcompoorga(f)+Mhumus(T1)+0.2*X[n]
    #X est la quantité d'engrais effective
    return Sol
#fdc le coefficient relatif à la gestion de la parcelle considérée
#Norganique représente le pourcentage d’azote total de la parcelle
#Da la densité apparente de l’horizon labouré 
# ep l’épaisseur de l’horizon labouré 
# %arg et CaCO3, respectivement les pourcentages d’argile et de calcaire de l’horizon labouré

#def CaU(mu,MS,MS7,SDT,SDT1,l) :
    #l et mu représentent respectivement l’ordonnée à l’origine et la pente de la relation linéaire entre le CAU et la vitesse de croissance de la culture dans les sept jours précédents l’apport [Ln(MSj/10) – ln(MSj-7/10)/(SDTj – SDTj-7)].
 #   CaU=min(100,l+(mu*100*(((ln(MS/10) – ln(MS7/10))*(MS/10))/(SDT–SDT7))))
  #  return CaU
T=[1,23,10,23,11,4,5,13,24,11,2,24,13,19,16,13,5,3,8,3]
X=[5,0,3,0,2,2,1,0,2,3,0,1,4,2,5,2,6,0,1,3]
print(sol(7,T,X))