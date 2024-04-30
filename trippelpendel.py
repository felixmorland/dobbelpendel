import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# True => animer bevegelsen
#anim = sys.argv[4]=='True'
anim = True

# Massene og lengdene til de tre pendelene
m1, m2, m3 = 1,1,1
L1, L2, L3 = 1,1,1

# Initialbetingelsene
#phi01, phi02, phi03 = [float(p) for p in sys.argv[1:4]]
phi01, phi02, phi03 = np.pi/2, 0, 0
dphi01, dphi02, dphi03 = 0,0,0

# Tyngdeakselerasjon
g = 9.81

# Tid
dt = 0.01
T = 10
t = np.arange(0,T,dt)
N = int(T/dt)

# Finne analytisk løsning for små utslag
omega = [0.41577, 2.2943, 6.2899]
Avkt = np.array([[(omega[i]**2-4*omega[i]+2)/(2*omega[i]) for i in range(3)],
                    [0.5*(2-omega[i]) for i in range(3)],
                    [1,1,1]])
C = np.linalg.solve(Avkt, np.array([phi01, phi02, phi03]))

anlphi = []
for i in range(3):
    phii = np.zeros(N)
    for j in range(3):
        phii += C[j]*Avkt[i,j]*np.cos(np.sqrt(omega[j]*(g/L1))*t)
    anlphi.append(phii)

# Analytiske løsninger
anlphi1, anlphi2, anlphi3 = anlphi
##################################

phi = np.zeros((N,3))   # Vinklene
dphi = np.zeros((N,3))  # Vinkelhastighetene
d2phi = np.zeros((N,3)) # Vinkelakselerasjon

phi[0] = phi01, phi02, phi03
dphi[0] = dphi01, dphi02, dphi03


for i in range(N-1):

    phi1, phi2, phi3 = phi[i,0], phi[i,1], phi[i,2]
    dphi1, dphi2, dphi3 = dphi[i,0], dphi[i,1], dphi[i,2]

    # A er en 3x3-matrise
    A = np.array([[(m1+m2+m3)*L1**2, (m2+m3)*L1*L2*np.cos(phi2-phi1), m3*L1*L3*np.cos(phi3-phi1)],
                 [(m2+m3)*L1*L2*np.cos(phi2-phi1), (m2+m3)*L2**2, m3*L2*L3*np.cos(phi3-phi2)],
                 [m3*L1*L2*np.cos(phi3-phi1), m3*L2*L3*np.cos(phi3-phi2), m3*L3**2]])

    # b er en tredimensjonal søylevektor
    b = np.array([[-(m1+m2+m3)*g*L1*np.sin(phi1) + (m2+m3)*L1*L2*np.sin(phi2-phi1)*dphi2**2 + m3*L1*L3*np.sin(phi3-phi1)*dphi3**2],
                 [-(m2+m3)*g*L2*np.sin(phi2) - (m2+m3)*L1*L2*np.sin(phi2-phi1)*dphi1**2 + m3*L2*L3*np.sin(phi3-phi2)*dphi3**2],
                 [-m3*g*L3*np.sin(phi3) - m3*L1*L3*np.sin(phi3-phi1)*dphi1**2 - m3*L2*L3*np.sin(phi3-phi2)*dphi2**2]])

    # d2phi (vektor) = A^-1 . b
    matprod = np.linalg.solve(A,b)

    d2phi[i] = matprod[0,0], matprod[1,0], matprod[2,0]

    # Euler-Cromer
    dphi[i+1] = dphi[i] + d2phi[i]*dt
    phi[i+1] = phi[i] + dphi[i+1]*dt


phi1, phi2, phi3 = phi[:,0], phi[:,1], phi[:,2]

def writedata(filnavn):
    with open(f'data/{filnavn}', 'w') as outfile:
        for i in range(N):
            outfile.write(f'{phi1[i]} {phi2[i]} {phi3[i]}\n')


# Posisjonene til de tre pendelene i (x,y)-planet
r1 = np.stack((L1*np.sin(phi1), -L1*np.cos(phi1)), axis=1)
r2 = r1 + np.stack((L2*np.sin(phi2), -L2*np.cos(phi2)), axis=1)
r3 = r2 + np.stack((L3*np.sin(phi3), -L3*np.cos(phi3)), axis=1)


# E = (
# 0.5*(m1+m2+m3)*(L1**2)*(dphi[:,0]**2) + 0.5*(m2+m3)*(L2**2)*(dphi[:,1]**2) + 0.5*m3*(L3**2)*(dphi[:,2]**2)
# + (m2+m3)*L1*L2*dphi[:,0]*dphi[:,1]*np.cos(phi[:,1]-phi[:,0]) + m3*L1*L3*dphi[:,0]*dphi[:,2]*np.cos(phi[:,2]-phi[:,0]) + m3*L2*L3*dphi[:,1]*dphi[:,2]*np.cos(phi[:,2]-phi[:,1])
# + (m1+m2+m3)*g*L1*(1-np.cos(phi[:,0])) + (m2+m3)*g*L2*(1-np.cos(phi[:,1])) + m3*g*L3*(1-np.cos(phi[:,2]))
# )
# plt.plot(t, E)
# plt.show()

plt.suptitle(f'phi0 (vektor) = [{phi01}, {phi02}, {phi03}]')

# Banen de tre pendelene tar
path1 = plt.plot(r1[0,0], r1[0,1], color='#ffcccc', markersize=0.1)
path2 = plt.plot(r2[0,0], r2[0,1], color='#ccffcc', markersize=0.1)
path3 = plt.plot(r3[0,0], r3[0,1], color='#ccccff', markersize=0.1)

# Stengene mellom massene
stang01 = plt.plot(np.array([0, r1[0,0]]), np.array([0, r1[0,1]]), 'k-')
stang12 = plt.plot(np.array([r2[0,0], r1[0,0]]), np.array([r2[0,1], r1[0,1]]), 'k-')
stang23 = plt.plot(np.array([r3[0,0], r2[0,0]]), np.array([r3[0,1], r2[0,1]]), 'k-')

# Posisjonen til punktmassene
pendel1 = plt.plot(r1[0,0], r1[0,1], color='#cd482c', marker='o', markersize=10*((3*m1)/(4*np.pi))**(1/3))
pendel2 = plt.plot(r2[0,0], r2[0,1], color='#42b349', marker='o', markersize=10*((3*m2)/(4*np.pi))**(1/3))
pendel3 = plt.plot(r3[0,0], r3[0,1], color='#427cb3', marker='o', markersize=10*((3*m3)/(4*np.pi))**(1/3))

# Definere aksene i (x,y)-planet
Lmax = L1+L2+L3
plt.xlim(-(Lmax+0.1),Lmax+0.1)
plt.ylim(-(Lmax+0.1),Lmax+0.1)
ax = plt.gca()
ax.set_aspect('equal', adjustable='box')
plt.xlabel('x [m]')
plt.ylabel('y [m]')

if anim:
    # Animasjonsløkke
    for i in range(N-1):

        # Oppdatere path
        path1[0].set_data(r1[:i,0], r1[:i,1])
        path2[0].set_data(r2[:i,0], r2[:i,1])
        path3[0].set_data(r3[:i,0], r3[:i,1])

        # Oppdatere stengene
        stang01[0].set_data(np.linspace(0, r1[i+1,0]), np.linspace(0, r1[i+1,1]))
        stang12[0].set_data(np.linspace(r2[i+1,0], r1[i+1,0]), np.linspace(r2[i+1,1], r1[i+1,1]))
        stang23[0].set_data(np.linspace(r3[i+1,0], r2[i+1,0]), np.linspace(r3[i+1,1], r2[i+1,1]))

        # Oppdatere posisjonen til punktmassene
        pendel1[0].set_data(r1[i+1,0], r1[i+1,1])
        pendel2[0].set_data(r2[i+1,0], r2[i+1,1])
        pendel3[0].set_data(r3[i+1,0], r3[i+1,1])

        plt.draw()
        # For å lage en .gif senere
        # if i % 3 == 0:
        #     plt.savefig(f'images/tmp_{i:05d}.png')
        plt.pause(dt)

#     os.system(
# """cd images
# convert -delay 1x100 tmp_*.png ../gifs/trippelpendel.gif
# rm tmp_*.png
# """
#     )
    plt.show()
