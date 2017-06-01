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
        IRRAD[k]=controlData["swv_dwn"][k]
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
    C=par['C']
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
    ebflo=par['ebflo']
    SENESCeb=par['SENESCeb']
    percNgrainflo=par['percNgrainflo']
    P1Gmax=par['P1Gmax']
    RDTmax=par['RDTmax']
    HP1Gmax=par['HP1Gmax']
    P1G0=par['P1G0']
    WSC=par['WSC']
    MSaa=par['MSaa']
    aremob=par['aremob']
    Rem1=par['Rem1']
    SENESC1=par['SENESC1']
    SENESC2=par['SENESC2']
    SENESC3=par['SENESC3']
    DC=par['DC']
    return(alpha,delta,beta,gamma,fsdc,Mr,Ma,mu,l,ebmax1,ebmax2,eimax,k,D,C,R,Vmax,Ncrit1,Ncrit2,Ncrit3,Nmax1,Nmax2,eb1,eb2,eb3,LAI1,LAI2,LAI3,RNG1,RNG2,NGmax,Tmin,Tmax,Toptimum,MSseuil,ebflo,SENESCeb,percNgrainflo,P1Gmax,RDTmax,HP1Gmax,P1G0,WSC,MSaa,aremob,Rem1,SENESC1,SENESC2,SENESC3,DC)

def attributionParametres2(exp):
    Qsolinit=exp['Qsolinit']
    flo = exp['flo']
    eWinter=exp['eWinter']
    dateepis1cm=exp['dateepis1cm']
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

    return(Qsolinit,flo,eWinter,dateepis1cm,Norg,Da,ep,pc_arg,CaCO3,Dapport1,Qapport1,Dapport2,Qapport2,Dapport3,Qapport3,Dapport4,Qapport4,Dapport5,Qapport5)

###############################################################################
##################            Initial state       #############################
###############################################################################

###STATE
#############
MS={0 : 178.3} #Le Moulon, Arche genotype, 2001, Matiere seche a la sortie de l'hiver
LAI={0:0}
SDJ={0:0}
SDT={0:0}
QN={0 : 0.1, 1 : 0.2}
Qsol={0:0}
INN_dicti={}

SUM_INN_Defficiency = 0
STemp = 0
Time_defficiency = 0
Global_Time_Defficiency=0
isDefficiencyYesterday = False
min_INN=1
Sf=0
precision=10**(-6)
NGM=0
N_n=QN[1]-QN[0]


MSflo=0
QNflo=0
SDTflo=0
QNGflo=0
LAIflo=0
QNaerflo=0
perc_Naer=0

S_JN=0
Hum={0:0}
MS1Gpot=dict()
MSGn=dict()
Vitrem=dict()
QNGn=dict()
QNaern=dict()
QNveg=dict()


###############################################################################
##################            SOL            ##################################
###############################################################################


def fmin(T,TSR):
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
    f=fmin(T2M[n],averageTSR())
    Sf+=f
    if n>=7:
        Sol = Mresidu(f,Sf) + Mcompoorga(f) + Mhumus(T2M[n]) + CaU(MS[n],MS[n-7],SDJ[n],SDJ[n-7])/100 * apport_engrais(n) - Q_plante_n+Qsol[n]
    else:
        if n!=0:
            Sol = Mresidu(f, Sf) + Mcompoorga(f) + Mhumus(T2M[n]) + Qsol[n] - Q_plante_n + 0.2*apport_engrais(n)
        else:
            Sol = Mresidu(f, Sf) + Mcompoorga(f) + Mhumus(T2M[n]) + Qsol[n]  +0.2*apport_engrais(n)

    # 0.2 est le CAU dans la literature
    Qsol[n+1] = Sol





###############################################################################
##################            PLANTE         ##################################
###############################################################################
def modelAzodyn(n):
    global MSflo
    global QNflo
    global SDTflo
    global QNGflo
    global LAIflo
    global QNaerflo
    if n<flo:
        res=calcul_avflor(n)
        if n==flo-1:
            MSflo=res
            MSGn[flo]=0
            QNflo=QN[flo]
            SDTflo=SDT[flo]
            QNGflo=percNgrainflo*P1Gm()/100
            QNGn[flo]=QNGflo
            LAIflo=LAI[flo]
            QNaerflo=QN[flo]/R
            QNaern[flo]=QNaerflo
            QNveg[flo]=QNaern[flo]-QNGn[flo]
    else:
        res=calcul_aprflo(n)
    return res
## AVANT FLORAISON
#module matière sèche
def calcul_avflor(n):
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
    global NGM
    MS1=MS.get(n)
    CalculSol(n,N_n)
    #test
    isnotSeuil = (MS1/1000) < MSseuil
    #On etablit la valeur des parametres de controle selon la date
    if n<=dateepis1cm:
        ebmax=ebmax1
    else:
        ebmax = ebmax2
    PAR=IRRAD[n]*C
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
    if n+1 == flo:
        NGM = nbr_grain(1 - min_INN, Global_Time_Defficiency)
        SDT[n+1]=JNMS(0,n+1)
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
def calcul_aprflo(n):

    global SDJ
    PAR = IRRAD[n] * C
    BMpotG(n,MS1Gpot,JNMS(int(flo),n))
     #on fait appel à la fonction du module avant floraison pour obtenir MSflo

    DJ = (T2MX[n] - T2MN[n]) / 2
    SDJ[n + 1] = DJ
    SDT[n+1]=JNMS(0,n+1)
    MS1 = f_cas5(n, MS[n], PAR,T2M[n])
    delta_MS=(MS1-MS[n])/1000
    QNinit=QNaer(n+1,delta_MS)
    QNveg[n+1]=QNinit-QNG(n+1)
    MSG(n)


    return ()


#Biomasse potentielle
def BMpotG(n,MS1Gpot,JN):
	#et à une humidité du sol à la capacité au champ), sera donné par un outil de calcul
    MS1Gpot[n]=(P1Gm()/1000)/(1+(((P1Gm()-P1G0)/P1G0)**((DC/2-JN)/(DC/2))))

#moisture function
def H(average_H):
    return 0.2+0.8*(average_H-0.1)/(0.7-0.6)

def Humidity_fake(n):
    Hum[n]=0.5
#Normalised day
def JN(n):
    # JN la durée en jours normalises de la phase floraison, (jours normalises : Un jour normalise correspond a un jour à une temperature de 15°C
    Humidity_fake(n)
    res=fmin(T2M[n],Toptimum)*H(Hum[n])
    return res

def JNMS(debut,fin):
    S=0
    for i in range(debut,fin+1):
        S+=JN(i)
    return S

#maximum mass of the grain
def P1Gm():
    return HP1Gmax*min(P1Gmax,1000*RDTmax/NGM)

#Accumulation réelle de matière sèche
def MSG(n): #f est le jour de floraison
    # QNG : dictionnaire contenant la quantité d’azote des grains le jour n
    # MSaa le coefficient d’estimation de la biomasse des grains issue du carbone transféré avec les acides aminés
    # MS.get(n) la matière sèche aérienne de la culture du jour n
    # MSflo la matière sèche aérienne de la culture à floraison, donc MS d'avant floraison le dernier jour
    # WSC0 la quantité de sucres solubles stockés à floraison et remobilisables
    # WSC le coefficient de remobilisation des sucres solubles.
    MSGn[n]=MSGn[n-1]+min(MS1Gpot[n]*NGM,((QNG[n]*MSaa)/10+(MS[n]-MSflo)/10 + (WSC0*WSC)/10))

def QNG(n): #il faut encore trouver les valeurs des constantes
    #Rem1 le coefficient de dégradation des protéines
    #JNMS j la somme des jours normalisés au jour j
    #La quantite d'azote remobilisable est fixe a la floraison
    QNrem=QNflo*aremob
    Vitrem[n]=Rem1*JNMS(0,n)
    if SDT[n]<(SDTflo+DC/2):
        QNGn[n]=QNGn[n-1]+min(MSGn[n-1]*0.028,QNrem*Vitrem.items()[1]+QNaern[n]-QNflo)
    else:
        QNGn[n]=QNGn[n-1]+QNrem*Vitrem+QNaern[n]-QNflo
    return QNGn[n]
#Accumulation de matiere seche et d'azote dans la culture entre floraison et recolte
def GLAI(n):
    '''LAIflo la surface foliaire de la culture à floraison
    SENESC liste des coefficients de réduction de la surface foliaire sous l’effet de la remobilisation de l’azote.
    QNveg représente la quantité d’azote des parties végétatives de la culture
    QNaerflo la quantité d’azote des parties aériennes de la culture à floraison'''
    LAI[n] = min(LAI[n-1],LAIflo*(SENESC1*math.log(SENESC2*(QNveg[n]/(QNaerflo-QNGflo))+SENESC3)))
    return LAI[n]

def QNaer(i,MS1):
    #DJ est la liste des degres-jours
    global perc_Naer
    sdj = 0
    for dj in (list(SDJ.items())[:i]):
        sdj += dj[1]
    if sdj < 200:
        if MS1/1000<MSseuil:
            QNaern[i] = QN[i-1]+ min(Qsol[i-1],min(MS1*Nmax1,Vmax*T2M[i]))
        else:
            QNaern[i] = QN[i - 1] + min(Qsol[i-1], min(MS1 * Nmax2, Vmax * T2M[i]))
    else:
        if perc_Naer==0:
            perc_Naer=QN[i-1]/MS1
        QNaern[i] = QN[i-1]+ min(Qsol[i-1],min(MS1*perc_Naer,Vmax*T2M[i]))
    return QNaern[i]


def f_cas5(n,x1,PAR,T):
    eb = calcul_eb_post_flo(n)
    LAIc = GLAI(n)
    y1 = 1 - math.exp(-k * LAIc)
    y2 = eb * eimax * PAR
    return x1 + f_croi(T) * y1 * y2


def calcul_eb_post_flo(n):
    if (QNveg[n]/QNaerflo)>=0.6:
        return ebflo
    elif(QNveg[n]/QNaerflo)>0.3:
        return max(ebflo*SENESCeb*(QNveg[n]/QNaerflo)-1,0)
    else:
        return 0
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
    alpha,delta,beta,gamma,fsdc,Mr,Ma,mu,l,ebmax1,ebmax2,eimax,k,D,C,R,Vmax,Ncrit1,Ncrit2,Ncrit3,Nmax1,Nmax2,eb1,eb2,eb3,LAI1,LAI2,LAI3,RNG1,RNG2,NGmax,Tmin,Tmax,Toptimum,MSseuil,ebflo,SENESCeb,percNgrainflo,P1Gmax,RDTmax,HP1Gmax,P1G0,WSC,MSaa,aremob,Rem1,SENESC1,SENESC2,SENESC3,DC=attributionParametres1(par)
    Qsol[0],flo, eWinter,dateepis1cm, Norg, Da, ep, pc_arg, CaCO3,Dapport1,Qapport1,Dapport2,Qapport2,Dapport3,Qapport3,Dapport4,Qapport4,Dapport5,Qapport5=attributionParametres2(exp)

    controlData = pd.read_csv(args.control_file, index_col=False, sep=",")
    DOY,IRRAD,T2M,T2MN,T2MX,RAIN,WS10M=lectureControle(controlData)

    #print(T2M)
    last_cycle = 310
    # print(args.parameters_file)
    for i in range(0,last_cycle-1):
        modelAzodyn(i)
        # print("time :", i, "| Matseche :", Matseche(i))
    DataSetMS = dictToList(MS)
    DataSetQsol = dictToList(Qsol)
    DataSet={}
    for i in range(0, len(DataSetMS)):
        DataSet[i] = str(DataSetMS[i][0])+",MS,"+str(DataSetMS[i][0]+1)+",MS,"+ str(DataSetMS[i][1])+",Qsol,"+str(DataSetQsol[i][1])
    #print(DataSet)
    df = pd.DataFrame(list(DataSet.items()), columns=['#id','observation_model,time,variable1,value1,variable2value2' ])
    df.to_csv(args.output_file,sep=",", index =False)#header=None, sep ="\t"
