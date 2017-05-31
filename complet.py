#!/usr/bin/env python3

# Il faut ajouter l'initialisation des dictionnaires
# Pour utiliser une dico (correction sur PARi) : pour acceder a un element soit dic[key] soit dic.get(n)
# la methode get permet de ne pas avoitr d'erreur si la cle n existe pas mais renvoi None


import argparse
import numpy as np
import pandas as pd
from scipy import *

###############################################################################
############################  hypothèse du programme   ########################
###############################################################################
'''
D: fraction interceptée supposé constante



'''

###############################################################################
##################            AUX            ##################################
###############################################################################

#fonction pour lecture des donnees des fichiers txt
def lectureDonnees(donnees):
    n = len(donnees)
    dictionnaire = {}
    for k in range(n):
        dictionnaire[donnees[0][k]] = donnees[1][k]
    return(dictionnaire)

#Fonction auxiliaire pour générer le fichier de sortie à partir d'un dictionnaire
def dictToList(dictionnaire):
    res = []
    for k in range(len(dictionnaire)):
        res.append((k,dictionnaire[k]))
    return(res)
def lectureControle(controlData):
    n = len(controlData)
    DOY = {} #day of the year
    IRRAD={} #Radiation
    T2M = {} #Average Air Temperature At 2 m Above The Surface Of The Earth (degreesC)
    T2MN={} #Min Ait Temperatire At 2 m Aboce The Surface Of The Earth(degreesC)
    T2MX = {}  # Max Ait Temperatire At 2 m Aboce The Surface Of The Earth(degreesC)
    RAIN = {} #Average Precipitation (mm/day)
    WS10M = {} #Wind Speed At 10 m Above The Surface Of The Earth (m/s)
    for k in range(n):
        DOY[k]=controlData["DOY"][k]-299
        IRRAD=controlData["swv_dwn"][k]
        T2M[k]=controlData["T2M"][k]
        T2MN[k] = controlData["T2MN"][k]
        T2MX[k] = controlData["T2MX"][k]
        RAIN[k]=controlData["RAIN"][k]
        WS10M[k]=controlData["WS10M"][k]
    return(DOY,IRRAD,T2M,T2MN,T2MX,RAIN,WS10M)


###PARAMETERS
#parameters = pd.read_csv(args.parameters_file)
#parameters = pd.read_csv(args.parameters_file, header = None, index_col=False, sep="\t")
#par = lectureDonnees(parameters)
def attributionParametres1(par):

    alpha=par["alpha"]
    delta=par['delta']
    beta=par['beta']
    gamma=par['gamma']
    fsdc=par['fsdc']
    Mr=par['Mr']
    Ma=par['Ma']
    mu=par['mu']
    l=par['l']
    ebmax1=par['ebmax1']
    ebmax2=par['ebmax2']
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
    eb1=par['eb1']
    eb2 = par['eb2']
    eb3 = par['eb3']
    LAI1 = par['LAI1']
    LAI2 = par['LAI2']
    LAI3 = par['LAI3']
    RNG1=par['RNG1']
    RNG2 = par['RNG2']
    NGmax=par['NGmax']
    Tmin=par['Tmin']
    Tmax=par['Tmax']
    Toptimum=par['Toptimum']
    MSseuil= par['MSseuil']

    return(alpha,delta,beta,gamma,fsdc,Mr,Ma,mu,l,ebmax1,ebmax2,eimax,k,D,R,Vmax,Ncrit1,Ncrit2,Ncrit3,Nmax1,Nmax2,eb1,eb2,eb3,LAI1,LAI2,LAI3,RNG1,RNG2,NGmax,Tmin,Tmax,Toptimum,MSseuil)

def attributionParametres2(exp):
    Qsolinit=exp['Qsolinit']
    flo = exp['flo']
    eWinter=exp['eWinter']
    Norg=exp['Norg']
    Da=exp['Da']
    ep=exp['ep']
    pc_arg=exp['arg']
    CaCO3=exp['CaCO3']
    Dapport1=exp['Dapport1']
    Qapport1=exp['Qapport1']
    Dapport2 = exp['Dapport2']
    Qapport2 = exp['Qapport2']
    Dapport3 = exp['Dapport3']
    Qapport3 = exp['Qapport3']
    Dapport4 = exp['Dapport4']
    Qapport4 = exp['Qapport4']
    Dapport5 = exp['Dapport5']
    Qapport5 = exp['Qapport5']

    return(Qsolinit,flo,eWinter,Norg,Da,ep,pc_arg,CaCO3,Dapport1,Qapport1,Dapport2,Qapport2,Dapport3,Qapport3,Dapport4,Qapport4,Dapport5,Qapport5)

###############################################################################
##################            Initial state       #############################
###############################################################################

###STATE
#############
MS={0 : 178.3} #Le Moulon, Arche genotype, 2001, Matiere seche a la sortie de l'hiver
LAI={0:0}
SDJ={0:0}
INN_dicti={}

SUM_INN_Defficiency = 0
STemp = 0
Time_defficiency = 0
Global_Time_Defficiency=0
isDefficiencyYesterday = False
min_INN=1
Sf=0
precision=10**(-6)
#############
# CORRECTION
#############
QN={0 : 0.1, 1 : 0.2, 2 : 0.3} #dict()
N_n=0.1
Qsol={0:0}


###CONTROL
#############
# CORRECTION
#############
Tmoyenne={1 : 7.5, 2 : 13.7} #value 1 corresponds to Oct-March whereas 2 for March-August
GR = {1 : 422,2 : 1602} #global radiation, periods divided as before
PARi = {1 : GR[1]*0.48 , 2: GR[2]*0.48 }


###############################################################################
##################            SOL            ##################################
###############################################################################


def fmin(n,T,TSR):
    # T, average temperature of day n
    #n, day number
    #TSR, Average temperature between the end of winter and crop
    if T < 0:
        return 0
    else:
        inter=math.exp(alpha*(T-TSR))
        return inter


def Mresidu(fmin,Sfmin):
    s = Mr * fmin/Sfmin
    return s

def Mcompoorga(fmin):
    s = Ma * fmin
    return s

def Mhumus(T):
    # pc arg représente le % d'argile
    if T < 0:
        return 0
    else:
        #print("Dans Mh")
        #print("fsdc",fsdc,";  expo:",math.exp(alpha * (T - Toptimum)), ";  reste",Norg/1000*Da*ep * (beta / ((pc_arg*gamma) * (CaCO3*delta))))
        Mh = fsdc * math.exp(alpha * (T - Toptimum)) /1000*Norg*Da*ep * (beta / ((pc_arg*gamma) * (CaCO3*delta)))
        return Mh

def CaU(MS,MS7,SDT,SDT7):
    #l et mu représentent respectivement l’ordonnée à l’origine et la pente de la relation linéaire entre le CAU et la vitesse de croissance de la culture dans les sept jours précédents l’apport [Ln(MSj/10) – ln(MSj-7/10)/(SDTj – SDTj-7)].
    CaU=min(100,l+(mu*(((math.log(MS/10)-math.log(MS7/10))*(MS/10))/(SDT-SDT7))))
    return CaU

def averageTSR():
    Somme=0
    long=len(T2M)
    long=int(long)
    for i in range(int(eWinter),long):
        Somme+=T2M[i]
    j=len(T2M)-int(eWinter)
    return Somme/j

def apport_engrais(n):
    if n==Dapport1:
        return Qapport1
    elif n==Dapport2:
        return Qapport2
    elif n==Dapport3:
        return Qapport3
    elif n==Dapport4:
        return Qapport4
    elif n==Dapport5:
        return Qapport5
    else:
        return 0

def CalculSol(n,Q_plante_n):
    global Sf
    f=fmin(n,T2M[n],averageTSR())
    Sf+=f
    #if n>=7:
    if n<0:
        print("CaU",CaU(MS[n],MS[n-7],SDJ[n],SDJ[n-7])/100)
        Sol = Mresidu(f,Sf) + Mcompoorga(f) + Mhumus(T2M[n]) + CaU(MS[n],MS[n-7],SDJ[n],SDJ[n-7])/100 * apport_engrais(n) - Q_plante_n+Qsol[n]
    else:
        if n!=0:
            Sol = Mresidu(f, Sf) + Mcompoorga(f) + Mhumus(T2M[n]) + Qsol[n] - Q_plante_n + 0.2*apport_engrais(n)
        else:
            Sol = Mresidu(f, Sf) + Mcompoorga(f) + Mhumus(T2M[n]) + Qsol[n]  +0.2*apport_engrais(n)
    '''
    print("n",n)
    print("Mr",Mresidu(f,Sf))
    print("Morg",Mcompoorga(f))
    print("Mh",Mhumus(T2M[n]))
    print("Qsol",Qsol[n])
    print("Qplante",Q_plante_n)
    print("Sol",Sol)
    print("  ")'''
    # X est la quantité d'engrais effective
    # P est la partie d'azote absorbée au temps j
    # 0.2 est le CAU dans la literature
    Qsol[n+1] = Sol





###############################################################################
##################            PLANTE         ##################################
###############################################################################
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
    global min_INN
    global Global_Time_Defficiency
    global SUM_INN_Defficiency
    global STemp
    global Time_defficiency
    global isDefficiencyYesterday
    global N_n
    MS1=MS.get(n)
    CalculSol(n,N_n)
    #test
    isnotSeuil = (MS1/1000) < MSseuil
    #On etablit la valeur des parametres de controle selon la date
    if n<141:
        PAR = PARi[1]
        Tmoy = Tmoyenne[1]
        ebmax=ebmax1
    else:
        PAR = PARi[2]
        Tmoy = Tmoyenne[2]
        ebmax = ebmax2

    #initialisation
    #Calcul INN
    isTneg=T2M[n]<0
    if isTneg:
        c1=0
    else:
        c1 = Vmax * T2M[n]

    if (n>1):
        delta_MS = (MS1 - MS[n-1])/1000
    else:
        delta_MS=0
    if isnotSeuil:
        #N_n_init = min(Qsol[n + 1], min(R * Nmax1 * (delta_MS), c1))
        N_n_init = min(Qsol[n + 1],c1)

        if N_n_init==0:
            INN=0
        else:
            INN = N_n_init / (delta_MS * Ncrit1)
    else:
        #N_n_init = min(Qsol[n + 1], min(R * Nmax2 * (delta_MS), c1))
        N_n_init = min(Qsol[n + 1],c1)
        if N_n_init==0:
            INN=0
        else:
            INN = N_n_init / (delta_MS * Ncrit2 * (MS1 / 1000) ** Ncrit3)
    N_n_loop=N_n_init
    #test
    DJ=(T2MX[n] - T2MN[n]) / 2
    SDJ[n + 1] = DJ
    #isDjok= (Dj>Tmin and DJ<Tmax)
    iscarence= INN<1
    whilecondition1=True
    whilecondition2= True
    arretwhile=True
    j=0
    while (whilecondition1 and whilecondition2 and arretwhile):
        #init
        N_n_init=N_n_loop
        #calcul
        #Sans carence
        if iscarence:
            Time_defficiency_loop = DJ+Time_defficiency
            SUM_INN_Defficiency_loop= INN+SUM_INN_Defficiency
            INNI_loop = SUM_INN_Defficiency_loop / Time_defficiency_loop
            if isnotSeuil:
                MS1,LAI_loop = f_cas3(n,MS1, PAR,ebmax,INNI_loop,T2M[n])
            # cas 4
            else:
                MS1, LAI_loop = f_cas4(n,MS1, PAR,ebmax,INNI_loop,T2M[n])
            isDefficiencyYesterday_loop=True
        else:
            #cas1
            if isnotSeuil:
                MS1, LAI_loop= f_cas1(MS1,PAR,ebmax,T2M[n])
            #cas 2
            else:
                MS1, LAI_loop=f_cas2(MS1,PAR,ebmax,T2M[n])
            isDefficiencyYesterday_loop=False

        #test
        delta_MS=(MS1-MS[n])/1000
        if isnotSeuil:
            N_n_loop = min(Qsol[n + 1], min(R * Nmax1 * (delta_MS), c1))

            if N_n_loop == 0:
                INN = 0
            else:
                INN = N_n_loop/ (delta_MS * Ncrit1)
        else:
            N_n_loop = min(Qsol[n + 1], min(R * Nmax2 * (delta_MS), c1))
            if N_n_loop == 0:
                INN = 0
            else:
                INN = N_n_loop/ (delta_MS * Ncrit2 * (MS1 / 1000) ** Ncrit3)
        if isnotSeuil:
            whilecondition1 = (N_n_init != min(Qsol[n + 1], min(R * Nmax1 * (delta_MS), c1)))
        else:
            whilecondition1 = (N_n_init != min(Qsol[n + 1], min(R * Nmax2 * (delta_MS), c1)))
        whilecondition2 = (N_n_init - N_n_loop > precision)
        j+=1
        arretwhile=j<10


    N_n=N_n_loop
    QN[n + 1] = QN[n] + N_n_loop
    LAI[n+1]=LAI_loop
    INN_dicti[n] = INN
    if isDefficiencyYesterday:
        Time_defficiency += DJ
        SUM_INN_Defficiency += INN
    else:
        Time_defficiency = 0
        SUM_INN_Defficiency = 0
    if (n <= flo):
        if INN < min_INN:
            min_INN = INN
        if INN < 0.9:
            Global_Time_Defficiency += DJ
    if n == flo:
        NGM = nbr_grain(1 - min_INN, Global_Time_Defficiency)
    MS[n+1]=MS1
    return MS1


def f_croi(T):
    if T<Tmin:
        return 0
    elif T>Toptimum:
        return 1-((T-Toptimum)/(Tmax-Toptimum))**2
    else:
        return 1 - ((T - Toptimum) / (Toptimum - Tmin)) ** 2


def f_cas1(x1,PAR,ebm,T):
    LAI=D*Ncrit1*x1
    y1=1-math.exp(-k*LAI)
    y2=ebm*eimax*PAR
    return x1+f_croi(T)*y1*y2,LAI

def f_cas2(x1,PAR,ebm,T):
    y=Ncrit2*(x1/1000)**Ncrit3
    LAI=D*y*x1
    y1=1-math.exp(-k*LAI)
    y2=ebm*eimax*PAR
    return x1+f_croi(T)*y1*y2,LAI

def f_cas3(n,x1,PAR,ebm,INNI,T):
    global LAI
    eb = calcul_eb_carence(ebm, INNI)
    LAIc=max(LAI[n],D*Ncrit1*x1*LAI1*(1-LAI2*math.exp(-LAI3*INNI)))
    y1=1-math.exp(-k*LAIc)
    y2=eb*eimax*PAR
    return x1+f_croi(T)*y1*y2,LAIc

def f_cas4(n,x1,PAR,ebm,INNI,T):
    global LAI
    eb=calcul_eb_carence(ebm,INNI)
    y=Ncrit2*x1**Ncrit3
    LAIc = min(LAI[n], D * y * x1 * LAI1 * (1 - LAI2 * math.exp(-LAI3 * INNI)))
    y1=1-math.exp(-k*LAIc)
    y2=eb*eimax*PAR
    return x1+f_croi(T)*y1*y2,LAIc

def calcul_eb_carence(ebma,INNI):
    return min(ebma,max(0,eb1*(1-eb2*math.exp(-eb3*INNI))))


def nbr_grain(IC, DC):
    RNG=min(1,RNG1-(RNG2*IC*DC))
    return (NGmax*NGmax*RNG)



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


    #QNGrains(3)
###############################################################################
##################            code           ##################################
###############################################################################
if __name__ == "__main__" :
    parser = argparse.ArgumentParser()
    parser.add_argument('--parameters-file', required=True)
    parser.add_argument('--experience-file', required=True)
    parser.add_argument('--control-file', required=True)
    parser.add_argument('--output-file', required=True)
    args = parser.parse_args()
    parameters = pd.read_csv(args.parameters_file, header = None, index_col=False, sep=",")
    par = lectureDonnees(parameters)
    experience = pd.read_csv(args.experience_file, header=None, index_col=False, sep=",")
    exp = lectureDonnees(experience)
    alpha,delta,beta,gamma,fsdc,Mr,Ma,mu,l,ebmax1,ebmax2,eimax,k,D,R,Vmax,Ncrit1,Ncrit2,Ncrit3,Nmax1,Nmax2,eb1,eb2,eb3,LAI1,LAI2,LAI3,RNG1,RNG2,NGmax,Tmin,Tmax,Toptimum,MSseuil=attributionParametres1(par)
    Qsol[0],flo, eWinter, Norg, Da, ep, pc_arg, CaCO3,Dapport1,Qapport1,Dapport2,Qapport2,Dapport3,Qapport3,Dapport4,Qapport4,Dapport5,Qapport5=attributionParametres2(exp)

    controlData = pd.read_csv(args.control_file, index_col=False, sep=",")
    DOY,IRRAD,T2M,T2MN,T2MX,RAIN,WS10M=lectureControle(controlData)
    #print(T2M)
    last_cycle = 200
    # print(args.parameters_file)
    for i in range(0,last_cycle-1):
        Matseche(i)
        # print("time :", i, "| Matseche :", Matseche(i))
    DataSetMS = dictToList(MS)
    DataSetQsol = dictToList(Qsol)
    DataSet=dict()
    for i in range(0, len(DataSetMS)):
        DataSet[i] = [DataSetMS[i][0], "MS", DataSetMS[i][0]+1, "MS", DataSetMS[i][1],"Qsol",DataSetQsol[i][1]]
    print(DataSet)
    df = pd.DataFrame(data = DataSet, columns=['#id','observation_model','time','variable1','value1','variable2','value2' ])
    df.to_csv(args.output_file, index =False)#header=None, sep ="\t"
