# -*- coding: utf-8 -*-
"""
Created on Tue Feb 14 15:59:28 2017
Lecture de fichiers
@author: elena
"""

import os
os.chdir("/Users/elena/Documents/2A/Projetinno/code") #ecrire le chemin où se trouve les fichiers qui vont être utilisés

def initialisationParametres(docname): #le nom du fichier doit inclure l'extension (.txt) et être écrit comme chaîne de caractères
    fichier = open(docname, "r")
    parametres = {}
    fichier.readline()
    for line in fichier.readlines():
        words = line.split()
        parametres[words[0]]=float(words[1]) #on utilise float pour accepter les paramètres à virgule
    fichier.close()
    return(parametres)
    
def initialisationVariablesControle(docname):
    fichier = open(docname, "r")
    variablescontrole = {}
    for line in fichier.readlines():
        words = line.split()
        variablescontrole[words[0]]=float(words[1]) #on utilise float pour accepter les paramètres à virgule
    fichier.close()
    return(variablescontrole)

#Pas vraiment utile d'en faire deux différents, mais comme c'est déjà codé...

#Une fois les données ont été estimées, cette fonction va permettre d'ouvrir le fichier docname, qui existe déjà, et copier les valeurs.
#Il pourrait être intéressant d'ajouter aussi une date ou un repère
def ecritureDonneesObtenues(donnees,docname): #les donnees sont sous forme de dictionnaire dont les cles sont des mots et les valeurs de float
    fichier = open(docname, "a") #a permet d'écrire à la fin et de ne pas écraser le fichier
    #si on désire créer un fichier, il faut changer "a" par "x"
    print("Quelle date donner à ce groupe de données?")    
    date = input("-->")
    fichier.write("\n"+date)
    for donnee in donnees.keys():
        valeur = donnees[donnee]        
        fichier.write("\n"+donnee+" "+str(valeur))
    fichier.close()
    
#Fonction pour la lecture de donnees à une date donnee en connaissant le nombre de données n à lire (pour arrêter la lecture et pas repartir sur d'autres dates)
def initialisationDonneesDatee(docname,date,n): #il faut faire attention à bien ecrire la date toujours sous le même format (jj mois aaaa)
    fichier = open(docname, "r")
    donnees = {}
    dateFichier = fichier.readline()
    while (dateFichier!=date): #on cherche la date sur le fichier
        dateFichier = fichier.readline()    
    for i in range(0,n): #on lit autant de lignes qu'il y a de donnees
        line = fichier.readline()
        words = line.split()
        donnees[words[0]]=float(words[1]) #on utilise float pour accepter les paramètres à virgule
    fichier.close()
    return(donnees)
    