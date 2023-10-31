# -*- coding: utf-8 -*-
#1D modell for krystall

import matplotlib.pyplot as plt
import numpy as np

#m=elektronmasse, V0 = potensialbrÃ¸nndybde
hbar=1.05E-34
m=9.11E-31
V0=-9.0*1.6E-19
#N = antall diskrete posisjoner i hver brÃ¸nn
N=20
#dz = skrittlengde
dz=1.0E-11
#Vi lager potensialet.
#Vunit er den repeterbare enheten med 1 brÃ¸nn og 1 barriere
Vunit = [V0]*N + [0]*N
#Nunit = antall brÃ¸nner (atomer)
Nunit = 50
V = [0.0]*4*N
for n in range(Nunit):
   V = V + Vunit
 
#Med 4*N til venstre og 3*N til hÃ¸yre blir hele systemet symmetrisk,
#med i alt Nunit brÃ¸nner (atomer)
V = V + [0.0]*3*N
#Ntot = antall elementer i V
Ntot = len(V)
#z = liste med posisjonsverdier (m)
z = [dz*n for n in range(Ntot)]
#d = liste med diagonalelementer i Hamiltonmatrisen H
d = [v + hbar**2/(m*dz**2) for v in V]
#e = verdi til ikke-diagonale elementer i H, dvs H(i,i+-1) = -e
e = - hbar**2/(2*m*dz**2)
#Initialisering av matrisen H: Legger inn verdi 0 i samtlige elementer
H = [[0]*(Ntot) for n in range(Ntot)]
#Dobbel for-lÃ¸kke som lager den tridiagonale matrisen H
for i in range(Ntot):
    for j in range(Ntot):
        if i==j:
            H[i][j]=d[i]
        if abs(i-j)==1:
            H[i][j]=e
#Finner w = egenverdiene og v = egenvektorene til matrisen H           
w,v = np.linalg.eigh(H)
#EineV = tabell med energiverdier i enheten eV
EineV = w/1.6E-19
#Neste 2 linjer skriver ut laveste energi og EineV[Nunit-1]
print('Laveste energi:',EineV[0])
print('Energi nr',Nunit,':',EineV[Nunit-1])
print('Laveste energi, bÃ¥nd 2:',EineV[Nunit])
print('Energi nr',2*Nunit,':',EineV[2*Nunit-1])
#VineV = tabell med potensialet i enheten eV
VineV = [x/1.6E-19 for x in V]
#Plotter |psi|^2 for tilstanden midt i energibÃ¥ndet
#plt.figure('Probability density')
#plt.plot(z,np.abs(v[:,int(Nunit/2-1)]**2))
#plt.title('1D crystal, probability density',fontsize=20)
#plt.xlabel('$z$ (m)',fontsize=20)
#plt.ylabel('$|\psi|^2$ for state in center of band',fontsize=20)
#plt.xlim(min(z),max(z))
#plt.show()
#Plotter V(z), samt de Nunit laveste energiverdiene som horisontale linjer
plt.figure('Energy band')
plt.plot(z,VineV)
for j in range(5*Nunit):
   l = plt.axhline(EineV[j])
   
plt.title('Energy band in 1D crystal',fontsize=20)
plt.xlabel('$z$ (m)',fontsize=20)
plt.ylabel('Potential $V$ and energies $E$ (eV)',fontsize=20)
plt.xlim(min(z),max(z))
plt.show()