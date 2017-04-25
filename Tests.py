#Je vais faire 2 tests, un avec un dico tout fait pour voir si ca tourne et l'autre avec la possibilité de le modifier
#première fonction sans dico à l'interieur ie on a les valeurs pour l'année à venir
dico
#celui-ci  c'est juste pour savoir si ca tourne
def test1(dico) :
    for i in range 365 :
        #365 c'est la durée de la recolte, a mieux estimer
        Tmoy=Tmoy(dico.T,i)
        Ms1=dico.Ms[i-1]+dico.ei+dico.eb
        Mn=dico.Mn+Mresidu(dico.Mr,fmin(i,T))+Mcompoorga(Ma.dico,fmin(i,T))+Mhumus(dico.fdsc,dico.a,dico.beta,dico.delta,Tmoy,dico.Tref,dico.Norga,dico.Da,dico.ep,dico.Beta,dico.pcarg,dico.CaC03,dico.d)+dico.Xi*CaU(dico.mu,dico.MS,dico.MS7,dico.SDT,dico.SDT1,dico.l)-dico.Pj
        dico.Mn=Mn
        #fin module sol
        
        
def test2(dico):
    for i in range 365:
        #on actualise le dico au tour par tour
        #il faut déterminer ce qui change CaCO3? pcarg?
        dico.compteur=i
        dico.T=dico.T+[input("Entrer la température au jour ",i)]
        dico.Xi=input("Entrer Xi")
        LAI=#Ca depend de comment on compte le calculer
        dico.ei=effintercept(dico.eimax,dico.k,LAI)
        dico.eb=econvlum(dico.ebmax,dico.Redeb)
        Ms1=dico.Ms[i-1]+dico.ei+dico.eb
        dico.Ms+=[Ms1]
        #SDT?
        Tmoy=Tmoy(dico.T,i)
        if i>7 :
            Ms7=dico.Ms[i-8]
        else :
            Ms7=0
        Mn1=Mn=dico.Mn+Mresidu(dico.Mr,fmin(i,T))+Mcompoorga(Ma.dico,fmin(i,T))+Mhumus(dico.fdsc,dico.a,dico.beta,dico.delta,Tmoy,dico.Tref,dico.Norga,dico.Da,dico.ep,dico.Beta,dico.pcarg,dico.CaC03,dico.d)+dico.Xi*CaU(dico.mu,Ms1,Ms7,dico.SDT,dico.SDT1,dico.l)-dico.Pj
        dico.Mn=Mn1
        #fin module sol
        
        