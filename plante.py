#!/usr/bin/env python3

# Il faut ajouter l'initialisation des dictionnaires
# Pour utiliser une dico (correction sur PARi) : pour acceder a un element soit dic[key] soit dic.get(n)
# la methode get permet de ne pas avoitr d'erreur si la cle n existe pas mais renvoi None


import argparse
import numpy as np
import pandas as pd
from scipy import *
from scipy.optimize import brentq

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
###PARAMETERS
#parameters = pd.read_csv(args.parameters_file)
#parameters = pd.read_csv(args.parameters_file, header = None, index_col=False, sep="\t")
#par = lectureDonnees(parameters)
def attributionParametres(par):
    ebmax=par['ebmax']
    eimax=par['eimax']
    k=par['k']
    D=par['D']
    R=par['R']
    Vmax=par['Vmax']
    #I have changed Ncrit into 3 coeff
    Ncrit1 = par['Ncrit1']
    Ncrit2 = par['Ncrit2']
    Ncrit3 = par['Ncrit3']
    Nmax1 = par['Nmax1']
    Nmax2 = par['Nmax2']
    MSseuil= par['MSseuil']
    flo = par['flo']
    return(ebmax,eimax,k,D,R,Vmax,Ncrit1,Ncrit2,Ncrit3,Nmax1,Nmax2,MSseuil,flo)

###STATE
#############
MS={0 : 178.3} #Le Moulon, Arche genotype, 2001, Matiere seche a la sortie de l'hiver
#############
# CORRECTION
#############
QN={0 : 0.1, 1 : 0.2, 2 : 0.3} #dict()


###CONTROL
#############
# CORRECTION
#############
Tmoyenne={1 : 7.5, 2 : 13.7} #value 1 corresponds to Oct-March whereas 2 for March-August
GR = {1 : 422,2 : 1602} #global radiation, periods divided as before
PARi = {1 : GR[1]*0.48 , 2: GR[2]*0.48 }
#############
# CORRECTION
#############
#LAI= dict()
#MS1Gpot = {}
#MSG={}
#QNG = {}
#QNaerien = {}
#DJ = list() #liste des degres-jours
#ei=list()
#eb=list()
#INNI=0
#INN0=0


#fonction pour lecture des donnees des fichiers txt
def lectureDonnees(donnees):
    n = len(donnees)
    dictionnaire = {}
    for k in range(n):
        dictionnaire[donnees[0][k]] = donnees[1][k]
    return(dictionnaire)
## AVANT FLORAISON

#module matière sèche
def Matseche(n):
    '''
    n: jour concerné que l'on va comparer à 141 pour savoir si on se trouve dans la periode oct-mars ou mars-aout
    MS: matiere séche dans la plante, dictionnaire contenant les jours av n
    PARi: rayonnement photosynthètique actif, depend de la periode concernee
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
    Nous pouvions étudier 4 cas de figure différents :
    Sans carence
        Cas 1 : MSj/1000<MSseuil
        Cas 2 : sinon
    Avec Carence
        Cas 3 : MSj/1000<MSseuil
        cas 4 : sinon
    '''
    #On etablit la valeur des parametres de controle selon la date
    if n<141:
        PAR = PARi[1]
        Tmoy = Tmoyenne[1]
    else:
        PAR = PARi[2]
        Tmoy = Tmoyenne[2]
    #HYPOTHESE 0:
    #On suppose qu'on a pas de carence azotée on boucle pour vérifier
    INN_0,INN_1=1,1

    #HYPOTHESE 1:
    #Obliger de remplacer MS(n) par MS(n-1) donc boucle while pour vérifier à
    #posteriori du calcul que l'hypothèse est valable

    #Initialisation de la boucle de vérification de l'hypothèse 0
    #initialisation de la boucle de vérification de l'hypothèse 1
    MS1,MS2=MS[n-1],MS[n-1]

    #Sans carence

    #cas1
    if (MS2/1000) < par['MSseuil']:
        MS1=MS2
        MS2= f_cas1(MS1,PAR)
        print(conv)
    #cas 2
    else:
        MS1=MS2
        MS2=f_cas2(MS1,PAR)
    #Boucle pour vérifier l'hypothèse 1
    #############
    # CORRECTION
    # Sseuil not define, je pense que vous vouliez mettre MSseuil
    #############
    # while ((MS1/1000<MSseuil and MS2/1000>=MSseuil) or (MS1/1000>=Sseuil and MS2/1000<MSseuil)):
    while ((MS1/1000<par['MSseuil'] and MS2/1000>=par['MSseuil']) or (MS1/1000>=par['MSseuil'] and MS2/1000<par['MSseuil'])):
        #cas1
        if (MS2/1000)<par['MSseuil']:
            f=lambda x: f_cas1(x,MS1,PAR)
            MS1=MS2
            MS2=brentq(f,MS1,MS1*2)
        #cas 2
        else:
            f=lambda x: f_cas2(x,MS1,PAR)
            MS1=MS2
            MS2=brentq(f,MS1,MS1*2)

    # Vérification de l'hypothèse 0 après l'initialisation
    delta_MS=MS2-MS[n-1]
    c1=Vmax*Tmoy
    #############
    # CORRECTION
    # Il manque la definition du dico Mn
    #############
    # Mn_n=Mn.get(n)
    # Mn_n=1
    # QN_n_1=QN.get(n-1)
    #############
    # CORRECTION
    # A on avis probleme d initialisation de MS du coup division par 0
    #############
    # if (MS2/1000)<MSseuil:
    #     N_n=min(Mn_n,min(R*Nmax[0]*(delta_MS),c1))
    #     QN[n]=QN_n_1+N_n
    #     INN_1=N_n/(delta_MS*Ncrit1)
    # else:
    #     N_n=min(Mn_n,min(R*Nmax[1]*(delta_MS),c1))
    #     QN[n]=QN_n_1+N_n
    #     INN_1=N_n/(delta_MS*Ncrit2*(MS2/1000)**Ncrit3)

    # Boucle de vérification de l'hypothèse 0
    while ((INN_0<1 and INN_1>=1) or (INN_0>=1 and INN_0<1)):
        if (INN_1<1):
            #Dernière boucle donne carence alors qu'avant non
            #cas3
            if (MS2/1000)<par['MSseuil']:
                f=lambda x: f_cas3(x,MS1,PAR)
                MS1=MS2
                MS2=brentq(f,MS1,MS1*2)
            #cas 4
            else:
                f=lambda x: f_cas4(x,MS1,PAR)
                MS1=MS2
                MS2=brentq(f,MS1,MS1*2)
            while (MS1/1000<par['MSseuil'] and MS2/1000>=par['MSseuil']) or (MS1/1000>=par['MSseuil'] and MS2/1000<par['MSseuil']):
                #cas3
                if (MS2/1000)<par['MSseuil']:
                    f=lambda x: f_cas3(x,MS1,PAR)
                    MS1=MS2
                    MS2=brentq(f,MS1,MS1*2)
                #cas 4
                else:
                    f=lambda x: f_cas4(x,MS1,PAR)
                    MS1=MS2
                    MS2=brentq(f,MS1,MS1*2)
        else:
            #Dernière boucle donne sans carence
            #cas1
            if (MS2/1000)<par['MSseuil']:
                f=lambda x: f_cas1(x,MS1,PAR)
                MS1=MS2
                MS2=brentq(f,MS1,MS1*2)
            #cas 2
            else:
                f=lambda x: f_cas2(x,MS1,PAR)
                MS1=MS2
                MS2=brentq(f,MS1,MS1*2)
            while (MS1/1000<par['MSseuil'] and MS2/1000>=par['MSseuil']) or (MS1/1000>=par['MSseuil'] and MS2/1000<par['MSseuil']):
                #cas1
                if (MS2/1000)<par['MSseuil']:
                    f=lambda x: f_cas1(x,MS1,PAR)
                    MS1=MS2
                    MS2=brentq(f,MS1,MS1*2)
                #cas 2
                else:
                    f=lambda x: f_cas2(x,MS1,PAR)
                    MS1=MS2
                    MS2=brentq(f,MS1,MS1*2)
        delta_MS=MS2-MS[n-1]
        INN_0=INN_1
        if (MS2/1000)<par['MSseuil']:
            N_n=min(Mn_n,min(R*Nmax1*(delta_MS),c1))
            QN[n]=QN_n_1+N_n
            INN_1=N_n/(delta_MS*Ncrit1)
        else:
            N_n=min(Mn_n,min(R*Nmax2*(delta_MS),c1))
            QN[n]=QN_n_1+N_n
            INN_1=N_n/(delta_MS*Ncrit2*(MS2/1000)**Ncrit3)
    #Avec Carence

    MS[n]=MS2
    return MS

def f_cas1(x1,PAR):
    y1=1-exp(-k*D*Ncrit1*x1)
    y2=ebmax*eimax*PAR
    return x1+y1*y2

def f_cas2(x1,PAR):
    y=Ncrit2*(x1/1000)**Ncrit3
    y1=1-exp(-k*D*y*x1)
    y2=ebmax*eimax*PAR
    return x1+y1*y2

def f_cas3(x,x1,PAR):
    y0=x-x1
    y1=1-exp(-k*D*Ncrit1*y0)
    y2=ebmax*eimax*PAR
    return y0-y1*y2

def f_cas4(x,x1,PAR):
    eb=calcul_eb_carence()
    y=Ncrit2*x**Ncrit3
    y0=x-x1
    y1=1-exp(-k*D*y*y0)
    y2=ebmax*eimax*PAR
    return y0-y1*y2

# def calcul_eb_carence():
#     Nperc=0

# def effintercept(eimax,k,LAI):
#     '''j'ai enlevé le LAI, on le calcule non ? NON LAI = D*QNcrit
#     '''
#     ei.append(eimax*(1-exp(-k*LAI)))
#     return

def Qazote(n, QN, MS, MN, R, Vmax, Tmoy, Nmax, MSeuil):
    #ici on va calculer la quantité d'azote accumulée dans la culture
    # MN est calculé par le module sol
    # Nmax contient les valeurs Nmax1 et Nmax2 qui donnent la teneur maximale en azote de la culture selon le niveau de biomasse du jour
    b = MS.get(n)/1000-MSeuil
    if b<0:
        perc_Nmax=Nmax1
    else :
        perc_Nmax=Nmax2
    QN[n]=QN.get(n-1)+ min(MN[n], min((R*MS[n-1]*perc_Nmax)-(MS[n]*perc_Nmax),Vmax*Tmoy))
    return()

# # Effets température et carence azote
# def econvlum(ebmax,Redeb):
#     '''
#     ebmax: constante permettant de calculer eb, on supposera eb=ebmax en 1ère approche
#             efficiene de conversion de la photosynthèse
#         =2.8 g.MJ−1
#     Redeb = liste des constantes necessaire pour le calcul en cas de carence
#         Redeb[1]= -3.13
#         Redeb[2]= 1.05
#         Redeb[3]= 2.18
#     '''

#     eb.append(min(ebmax,Redeb[0]*(1-Redeb[1]*exp(-Redeb[2]*INNI[n]))))
#     return()

# def INNI(Ncrit,bm_booleen):
#     '''
#     Cette fonction sert à caculer l'indice de nutrition azotée intégré INNI
#     Ncrit : liste de constante utile en cas de carrence
#         Ncrit1= 4.4
#         Ncrit2= 5.35
#         Ncrit3= -0.442
#     bm_booleen: booléen permettant d'ajuster le calcul
#     %N a la floraison = 2.8 (Girard, 1997)
#     '''
#     if bm_booleen==false:
#         perc_Ncrit=Ncrit[0]
#     else :
#         perc_Ncrit=Ncrit[1]*(MS(n-1)/1000)^(Ncrit[2])
#     perc_N=0 # mesuré
#     INN1=(perc_N/perc_Ncrit)
#     INNI=0
#     return INNI

# def LAIcarence(n,LAI,LAIc,RedLAI,INNI):
#     # n est le jour considéré
#     # LAI l'indice foliaire sans carence
#     # LAI l'indice foliaire avec carence
#     # RedLAI la liste des coefficients de reduction d'expansion foliaire
#     # INNI est l'indice de nutrition azotée intégrée
#     LAIc[n]=min(LAIc.get(n-1), LAI(n)*RedLAI[0]*(1-(RedLAI[1]*exp(-RedLAI[2]*INNI.get(n)))))
#     return

# def nombreDeGrains(RedNG,IC,DC):
#     #RedNG contient les coefficients
#     #IC représente l'intensité de la carence azotée
#     #DC la durée de la carence
#     RNG = min(1,RedNG[0]-(RedNG[1]*IC*DC))
#     return(RNG)

##APRES FLORAISON
#Calcul de l'azote accumulé dans les grains au jour n, WSC0
#Il manque : calcul JN, MSflo (quel est le jour de la floraison), QNflo, Rem1, JNMS
def QNGrains(n):
    '''
    n: jour concerné
    MSG: matière séche dans les grains, dictionnaire contenant les jours av n
    QNflo : quantité d'azote à la floraison, calculée grâce au module précédent
    aremob : coefficient à estimer grâce à littérature
    Rem1 le coefficient de dégradation des protéines
    JNMSn la somme des jours normalisés au jour n.
    eb : liste des efficiences de conversion de la culture (tableau où indice = jour)
    ebflo : efficience de conversion du rayonnement à floraison (obtenue grâce au module avant floraison)
    SENESCeb : coefficient de réduction de l’efficience de conversion du fait des transferts d’azote vers les grains.
    DC : durée maximale de la phase de remplissage des grains en jours normalisés
    QNcrit
    Nous pouvons étudier 2 cas de figure différents pour l'accumulation d'azote dans les grains :
        Cas 1 : remplissage des grains (SDTn < SDT flo + DC/2)
        Cas 2 : demande en azote des grains n'intervient plus
    Nous pouvons aussi étudier deux cas pour l'accumulation d'azote dans la culture selon la somme des degrés-jours depuis le semis au jour j :
    	Cas 1 : SDJn < 200
    	Cas 2 : SDJn >= 200
    	Pour un jour donné, on a DJ = (Tmin+Tmax)/2, la température de base est considérée nulle pour le blé
    '''
    #données dont on a besoin
    #jour de floraison :
    MSflo = Matseche(flo)[flo]
    BMpot(n,MS1Gpot,JN.get(n))
    MSG(n,MSG,MS1Gpot,QNG,MSaa,MSflo,WSC0) #on fait appel à la fonction du module avant floraison pour obtenir MSflo
    QNG(n,QNG,MSG,QNaerien,QNflo,Rem1,JNMS)

    return ()


# #Biomasse potentielle
# def BMpot(n,MS1Gpot,JN):
#     #P1G max le poids d’un grain maximum, littérature : P1G max(16% H) 48.5 (Gate, 1995)
#     P1Gmax = 48.5
#     #P1G 0 le poids d’un grain à floraison, littérature : (Girard, 1997)
#     P1G0 = 0.6
#     #DC représente la durée maximale de la phase de remplissage des grains en jours normalisés
#     DC = 90
#     #JN la durée en jours normalisés de la phase floraison, (jours normalisés : Un jour normalisé correspond à un jour à une température de 15°C
# 	#et à une humidité du sol à la capacité au champ), sera donné par un outil de calcul
#     MS1Gpot[n]=(P1Gmax/1000)/(1+(((P1Gmax-P1G0)/P1G0)**((DC/2-JN.get(n))/(DC/2))))
#     return()

# #Accumulation réelle de matière sèche
# def MSG(n,MSG,MS1Gpot,QNG,MSaa,MS,f,WSC0): #f est le jour de floraison
#     # QNG : dictionnaire contenant la quantité d’azote des grains le jour n
#     # MSaa le coefficient d’estimation de la biomasse des grains issue du carbone transféré avec les acides aminés
#     # MS.get(n) la matière sèche aérienne de la culture du jour n
#     # MSflo la matière sèche aérienne de la culture à floraison, donc MS d'avant floraison le dernier jour
#     MSflo = MS.get(f)
#     # NGM2 : est calcul par une fonction auxiliaire
#     NGM2 = NGM2()
#     # WSC0 la quantité de sucres solubles stockés à floraison et remobilisables
#     # WSC le coefficient de remobilisation des sucres solubles.
#     WSC = 0.75
#     QNG(n,QNG,MSG,QNaerien,QNflo,d)
#     MSG[n]=MSG.get(n-1)+min(MS1Gpot.get(n)*NGM2,((QNG.get(n)*MSaa)/10+(MS.get(n)-MSflo)/10 + (WSC0*WSC)/10))
#     return()

# def QNG(n,QNG,MSG,QNaerien,QNflo,Rem1,JNMS): #il faut encore trouver les valeurs des constantes
#     #Rem1 le coefficient de dégradation des protéines
#     #JNMS j la somme des jours normalisés au jour j
#     #La quantite d'azote remobilisable est fixe a la floraison
#     aremob = 0.68 #donnée littérature
#     QNrem=QNflo*aremob
#     Vitrem[n]=Rem1*JNMS.get(n)
#     if SDT.get(n)<(SDTflo+DC/2):
#         QNG[n]=QNG.get(n-1)+min(MSG.get(n)*0.028,QNrem*Vitrem+(QNaerien.get(n)-QNflo))
#     else:
#         QNG[n]=QNG.get(n-1)+QNrem*Vitrem+(QNaerien.get(n)-QNflo)
#     return()

# def NGM2():
# 	#on va utiliser les données de la littérature pour RedNG1 et RedNG2:
# 	RedNG1 = 1.00355
# 	RedNG2 = 0.0011
# 	IC = 0
# 	DC = 90
# 	NGmax = 27000
# 	RNG =min(1,RedNG1- (RedNG2*IC*DC))
# 	return(NGmax*NGmax*RNG)

# #Accumulation de matiere seche et d'azote dans la culture entre floraison et recolte
# #def GLAI(n,LAIflo,SENESC,QNveg,QNaerflo,QNGflo):
#     # LAIflo la surface foliaire de la culture à floraison
#     # SENESC liste des coefficients de réduction de la surface foliaire sous l’effet de la remobilisation de l’azote.
#     # QNveg représente la quantité d’azote des parties végétatives de la culture
#     # QNaerflo la quantité d’azote des parties aériennes de la culture à floraison
#  #   GLAI[n] = min(GLAI.get(n−1), LAIflo * (SENESC[0] * ln(SENESC[1] * (QNveg.get(n) /(QNarflo − QNGflo)) + SENESC[2])))

# def QNaerien(n,QNaerien,QN,MS):
# 	#DJ est la liste des degres-jours
# 	sdj = 0
# 	for dj in DJ[:n]:
# 		sdj += dj
# 	if sdj < 200 :
# 		QNaerien[n] = QN.get(n-1)+ min(MS.get(j)*pour_Nmax,0.5*Tmoy[n])
# 	else:
# 		QNaerien[n] = QN.get(n-1)+ min(MS.get(j)*pour_N[n-1]*Tmoy[n])
# 	return()

#Fonction auxiliaire pour générer le fichier de sortie à partir d'un dictionnaire
def dictToList(dictionnaire):
    res = []
    for k in range(len(dictionnaire)):
        res.append((k,dictionnaire[k]))
    return(res)
def lectureControle(controlData):
    n = len(controlData)
    DOY = {} #day of the year
    T2M = {} #Average Air Temperature At 2 m Above The Surface Of The Earth (degreesC)
    RAIN = {} #Average Precipitation (mm/day)  
    WS10M = {} #Wind Speed At 10 m Above The Surface Of The Earth (m/s)  
    for k in range(n):
        DOY[k]=controlData["DOY"][k]-299
        T2M[k]=controlData["T2M"][k]
        RAIN[k]=controlData["RAIN"][k]
        WS10M[k]=controlData["WS10M"][k]
    return(DOY,T2M,RAIN,WS10M)

if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--parameters-file', required=True)
    parser.add_argument('--control-file', required=True)
    parser.add_argument('--output-file', required=True)
    args = parser.parse_args()
    parameters = pd.read_csv(args.parameters_file, header = None, index_col=False, sep="\t")
    par = lectureDonnees(parameters)
    ebmax,eimax,k,D,R,Vmax,Ncrit1,Ncrit2,Ncrit3,Nmax1,Nmax2,MSseuil,flo=attributionParametres(par)
    controlData = pd.read_csv(args.control_file, index_col=False, sep=",")
    DOY,T2M,RAIN,WS10M=lectureControle(controlData)
    #print(T2M)
    last_cycle = 20
    print(args.parameters_file)
    for i in range(1,last_cycle):
        print("time :", i, "| Matseche :", Matseche(i))
    #DataSet = dictToList(Matseche(last_cycle))
    #df = pd.DataFrame(data = DataSet, columns=['id','observation_model','time','variable','value' ])
    #df.to_csv(args.output_file, index =False)#header=None, sep ="\t" 

    #QNGrains(3)
