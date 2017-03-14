#module matière sèche 
def Matseche(n,MS,eimax,k,LAI,ebmax,D,QNcrit):
    ei=effintercept(eimax,k,LAI,D,QNcrit)
    eb=econvlum(ebmax,k,LAI,D,QNcrit)
    MS+=[MS(n-1)+ei[n]+eb[n]]
    return MS
def effintercept(eimax,k,LAI,D,QNcrit):
    