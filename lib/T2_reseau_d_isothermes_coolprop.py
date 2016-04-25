#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Sauf mention explicite du contraire par la suite, ce travail a été fait par 
# Jean-Julien Fleck, professeur de physique/IPT en PCSI1 au lycée Kléber. 
# Vous êtes libres de le réutiliser et de le modifier selon vos besoins.

import numpy as np               # Les outils mathématiques
import CoolProp.CoolProp as CP   # Les outils thermodynamiques
import matplotlib.pyplot as plt  # Les outils graphiques
import os
os.chdir("/home/wilk/git-repositories/py4phys/lib/")
def isothermes_d_andrews(fluide,dico={}):
    """ Dessines les isothermes d'Andrews pour le fluide demandé avec des 
    choix par défaut qui peuvent être "overridden" en spécifiant ceux à 
    changer dans le dictionnaire 'dico'. Les options disponibles sont:
    * 'vmin' et 'vmax' pour définir les limites des échantillonnages en volume 
    massique. Par défaut vtripleL et vtripleG * 10
    * 'Prange' pour l'intervalle de pression affiché
    * 'T': une liste des températures pour lesquelles il faut tracer 
    l'isotherme.
    * 'titre': le titre à donner au graphique.
    * 'fichier': le nom du fichier dans lequel enregistrer la figure.
    * 'logx': Booléen indiquant si on veut un axe logarithmique en abscisse
    * 'logy': Booléen indiquant si on veut un axe logarithmique en ordonnée
    * 'legend': Booléen indiquant si on veut rajouter les légendes
    * 'saturation': Booléen indiquant si on veut rajouter la courbe de 
    saturation au tracé (défaut à False)
    """
    Pcritique = .9999*CP.PropsSI(fluide,'pcrit')  # Pression
    Tcritique = .9999*CP.PropsSI(fluide,'Tcrit')  # et température critique
    # Ptriple = 1.001*CP.PropsSI(fluide,'ptriple')  # Pression 
    # Ttriple = 1.001*CP.PropsSI(fluide,'Ttriple')  # et température au point triple
    Ptriple = CP.PropsSI(fluide,'ptriple')  # Pression 
    Ttriple = CP.PropsSI(fluide,'Ttriple')  # et température au point triple
    # On récupère les volumes massiques via les 'densités' (ie masses 
    # volumiques) données par CoolProp
    vtripleL = 1/CP.PropsSI('D','P',Ptriple,'Q',0,fluide)
    vtripleG = 1/CP.PropsSI('D','P',Ptriple,'Q',1,fluide)
    vcritique= 1/CP.PropsSI('D','P',Pcritique,'T',Tcritique,fluide)
    # L'ensemble des valeurs par défaut.
    DEFAUTS = {'vmin':vtripleL, 'vmax':vtripleG*10, 
       'Prange': None, 
       'T': np.arange(Ttriple,Tcritique*1.2,20),
       'titre': "Isotherme d'Andrews pour le fluide {}".format(fluide),
       'fichier': 'PNG/T2_reseau_d_isothermes_coolprop_{}.png'.format(fluide),
       'logx': True, 'logy': True, 'legend': True,
       'export-points': False,
       'fichier-Pv': 'PNG/T2_reseau_d_isothermes_coolprop_Pv{}.png'.format(fluide),
       'fichier-sat': 'PNG/T2_reseau_d_isothermes_coolprop_sat{}.png'.format(fluide),
       'saturation': False}
    DEFAUTS.update(dico)      # Mise à jour des valeurs par défaut via 'dico' 
    # L'échantillonnage sera différent
    if DEFAUTS['logx']:       # si l'axe est logarithmique
       v=np.logspace(np.log10(DEFAUTS['vmin']),np.log10(DEFAUTS['vmax']),600)
    else:                     # ou simplement linéaire
       v=np.linspace(DEFAUTS['vmin'],DEFAUTS['vmax'],200)
    if DEFAUTS['export-points']:
        f=open(DEFAUTS['fichier-Pv'],'ab') # Pour vider le fichier
        os.ftruncate(f.fileno(),0)
        os.lseek(f.fileno(),0, os.SEEK_SET)
    for Ti in DEFAUTS['T']:   # Tracé des différentes isothermes
        P = CP.PropsSI('P','T',Ti,'D',1/v,fluide)
        if DEFAUTS['export-points']:
            # np.savetxt(f,[[v[i],P[i]] for i in range(len(v))],delimiter='\t',header='v(m^3/kg),P(Pa)')
            np.savetxt(f,list(zip(v,P/1e5)),delimiter='\t',header='v(m^3/kg),P(bar)')
            np.savetxt(f,[np.nan],delimiter='\t')
        plt.plot(v,P,label='$T={}$'.format(Ti))
    if DEFAUTS['export-points']:
        f.close()
    if DEFAUTS['saturation']: # Tracé de la courbe de saturation
        P_sat= np.linspace(Ptriple,Pcritique,600)
        v_eb   = 1/CP.PropsSI('D','P',P_sat,'Q',0,fluide)
        v_rosee= 1/CP.PropsSI('D','P',P_sat,'Q',1,fluide)
        if DEFAUTS['export-points']:
            f=open(DEFAUTS['fichier-sat'],'ab')
            os.ftruncate(f.fileno(), 0)
            os.lseek(f.fileno(),0, os.SEEK_SET)
            # np.savetxt(f,[[P_sat[i],v_eb[i],v_rosee[i]] for i in range(len(P_sat))],header='P(Pa),v_eb(m^3/kg),v_rosee(m^3/kg)')
            np.savetxt(f,list(zip(P_sat/1e5,v_eb,v_rosee)),delimiter='\t',header='P(bar),Veb,Vrosee(m^3/kg)')
            f.close()
        plt.plot(v_eb,P_sat,'k',linewidth=2.0)
        plt.plot(v_rosee,P_sat,'k',linewidth=2.0)
    if DEFAUTS['Prange']: plt.ylim(DEFAUTS['Prange']) # Intervalle vertical
    plt.xlim((DEFAUTS['vmin'],DEFAUTS['vmax']))       # Intervalle horizontal
    if DEFAUTS['logx']: plt.xscale('log')             # Echelle log en x
    if DEFAUTS['logy']: plt.yscale('log')             # Echelle log en y
    if DEFAUTS['legend']: plt.legend()                # Rajout des légendes
    plt.xlabel('Volume massique $v$ en m$^3/$kg')     # Légende en abscisse
    plt.ylabel('Pression en Pa')                      # Légende en ordonnée
    plt.title(DEFAUTS['titre'])                       # Titre
    plt.savefig(DEFAUTS['fichier'])                   # Enregistrement
    plt.clf()                                         # Nettoyage

# Le fluide à étudier (à choisir parmi ceux donnés par CP.FluidsList())
fluide = 'carbondioxide'

# Le diagramme "par défaut"
isothermes_d_andrews(fluide)

# Les valeurs suivantes ont été choisies suite à l'observation du diagramme 
# par défaut. Il faudra certainement changer les valeurs si vous modifiez le 
# fluide
dico = {'Prange':(1e6,1e7),
        'fichier':'PNG/T2_reseau_d_isothermes_coolprop_{}_lin.png'.format(fluide),
        'logx':False, 'logy': False,
        'vmin': 1e-3, 'vmax':1e-2,
        # 'T': [200 + i*5 for i in range(30)],
        'T': [270,290,304.2,310,320], 
        'legend': False,
        'export-points': True,
        'fichier-Pv': '/home/wilk/enseignement/schemas/thermodynamique/coolprop/isothermes-andrews-Pv-{}.csv'.format(fluide),
       'fichier-sat': '/home/wilk/enseignement/schemas/thermodynamique/coolprop/isothermes-andrews-sat-{}.csv'.format(fluide),
        }
isothermes_d_andrews(fluide,dico)

# Le même en rajoutant la courbe de saturation
dico['saturation'] = True
dico['fichier'] = 'PNG/T2_reseau_d_isothermes_coolprop_{}_lin_sat.png'.format(fluide)
isothermes_d_andrews(fluide,dico)



