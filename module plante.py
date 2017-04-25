from math import exp

###############################################################################
############################  hypothèse du programme   ########################
###############################################################################
'''
D: fraction interceptée supposé constante



'''
###############################################################################
##################            code           ##################################
###############################################################################
#initialisation

ei=list()
eb=list()
INNI=0
INN0=0

## AVANT FLORAISON

#module matière sèche 
def Matseche(n,MS,PARi,eimax,k,ebmax,LAI):
    '''j'ai enlevé le LAI, on le calcule non ? NON
    n: jour concerné
    MS: matière séche dans la plante, dictionnaire contenant les jours av n - clés = entier, valeurs = chiffres
    PARi: rayonnement photosynthètique actif, dictionnaire contenant les n jours
        PARi[n]=D*RG[n]
    RG :  rayonnement global incident, dictionnaire
    D: fraction interceptée supposé constante
    eimax : constante permettant de calculer ei, inutile puisque constante
            ei étant l'efficacité d'interception de la culture
        =0.96 (Jeuffroy et Recous, 1999)
    k : paramètres variant selon l'espèce et la culture (de l'angle)
        =radiation extinction coefficient
        =0.72 (Jeuffroy et Recous, 1999)
    ebmax: constante permettant de calculer eb, on supposera eb=ebmax en 1ère approche
            efficiene de conversion de la photosynthèse
        =2.8 g.MJ−1
    D: fraction reliant LAI à QNcr
        =0.028 (valeur venant de biblio inconnue)
    QNcrit
    '''
    effintercept(eimax,k,LAI)
    econvlum(ebmax,k,LAI)
    MS[n]=[MS.get(n-1)+ei[n]+eb[n]*PARi[n]]  #ça serait plutôt +eb*PARi et on crée la nouvelle valeur dans MS associee a n
    return MS.get(n)
def effintercept(eimax,k,LAI): 
    '''j'ai enlevé le LAI, on le calcule non ? NON LAI = D*QNcrit
    '''
    ei.append(eimax*(1-exp(-k*LAI)))
    return

def Qazote(n, QN, MS, MN, R, Vmax, Tmoy, Nmax, MSeuil):
    #ici on va calculer la quantité d'azote accumulée dans la culture
    # MN est calculé par le module sol
    # Nmax contient les valeurs Nmax1 et Nmax2 qui donnent la teneur maximale en azote de la culture selon le niveau de biomasse du jour
    b = MS.get(n)/1000-MSeuil 
    if b<0:
        perc_Nmax=Nmax[0]
    else :
        perc_Nmax=Nmax[1]
    QN[n]=QN.get(n-1)+ min(MN.get(n), min((R*MS.get(n-1)*perc_Nmax)-(MS.get(n)*perc_Nmax),Vmax*Tmoy))
    return()

# Effets température et carence azote
def econvlum(ebmax,Redeb):
    '''
    ebmax: constante permettant de calculer eb, on supposera eb=ebmax en 1ère approche
            efficiene de conversion de la photosynthèse
        =2.8 g.MJ−1
    Redeb = liste des constantes necessaire pour le calcul en cas de carence
        Redeb[1]= -3.13
        Redeb[2]= 1.05
        Redeb[3]= 2.18
    '''
    
    eb.append(min(ebmax,Redeb[0]*(1-Redeb[1]*exp(-Redeb[2]*INNI[n]))))
    return()

def INNI(Ncrit,bm_booleen):
    '''
    Cette fonction sert à caculer l'indice de nutrition azotée intégré INNI
    Ncrit : liste de constante utile en cas de carrence
        Ncrit1= 4.4
        Ncrit2= 5.35
        Ncrit3= -0.442
    bm_booleen: booléen permettant d'ajuster le calcul
    %N a la floraison = 2.8 (Girard, 1997)
    '''
    if bm_booleen==false:
        perc_Ncrit=Ncrit[0]
    else :
        perc_Ncrit=Ncrit[1]*(MS(n-1)/1000)^(Ncrit[2])
    perc_N= # mesuré
    INN1=(perc_N/perc_Ncrit)
    INNI=
    return INNI

def LAIcarence(n,LAI,LAIc,RedLAI,INNI):
    # n est le jour considéré
    # LAI l'indice foliaire sans carence
    # LAI l'indice foliaire avec carence
    # RedLAI la liste des coefficients de reduction d'expansion foliaire
    # INNI est l'indice de nutrition azotée intégrée
    LAIc[n]=min(LAIc.get(n-1), LAI(n)*RedLAI[0]*(1-(RedLAI[1]*exp(-RedLAI[2]*INNI.get(n)))))
    return

def nombreDeGrains(RedNG,IC,DC):
    #RedNG contient les coefficients 
    #IC représente l'intensité de la carence azotée
    #DC la durée de la carence
    RNG = min(1,RedNG[0]-(RedNG[1]*IC*DC))
    return(RNG)
    
    ##APRES FLORAISON
    
    
