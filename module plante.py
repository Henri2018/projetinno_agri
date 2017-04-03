#module matière sèche 
def Matseche(n,MS,eimax,k,ebmax,D,QNcrit):#j'ai enlevé le LAI, on le calcule non ?
    ei=effintercept(eimax,k,D,QNcrit)
    eb=econvlum(ebmax,k,D,QNcrit)
    MS+=[MS(n-1)+ei[n]+eb[n]]
    return MS
def effintercept(eimax,k,D,QNcrit): #j'ai enlevé le LAI, on le calcule non ?
    LAI = D*QNcrit
    ei = eimax*(1-exp(-k*LAI))
    return ei
def econvlum(ebmax,k,D,QNcrit):
    
