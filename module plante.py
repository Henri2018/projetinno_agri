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
#module matière sèche 
def Matseche(n,MS,PARi,eimax,k,ebmax,LAI):
    '''j'ai enlevé le LAI, on le calcule non ? NON
    n: jour concerné
    MS: matière séche dans la plante, dictionnaire contenant les jours av n
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
    MS+=[MS(n-1)+ei[n]*eb[n]*PARi[n]]
    return MS
def effintercept(eimax,k,LAI): 
    '''j'ai enlevé le LAI, on le calcule non ? NON LAI = D*QNcrit
    '''
    ei.append(eimax*(1-exp(-k*LAI)))

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
    
