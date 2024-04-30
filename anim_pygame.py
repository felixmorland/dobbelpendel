import numpy as np
import pygame
import sys

# Massene og lengdene til de tre pendelene
m1, m2, m3 = 1, 1, 1
L1, L2, L3 = 1, 1, 1

# Initialbetingelsene
phi01, phi02, phi03 = np.pi / 2, 0, 0
dphi01, dphi02, dphi03 = 0, 0, 0

# Tyngdeakselerasjon
g = -9.81

# Tid
dt = 0.01
T = 10
t = np.arange(0, T, dt)
N = int(T / dt)

# Posisjonene til de tre pendelene i (x, y)-planet
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

r1 = np.stack((L1*np.sin(phi[:,0]), -L1*np.cos(phi[:,0])), axis=1)
r2 = r1 + np.stack((L2*np.sin(phi[:,1]), -L2*np.cos(phi[:,1])), axis=1)
r3 = r2 + np.stack((L3*np.sin(phi[:,2]), -L3*np.cos(phi[:,2])), axis=1)

def calculate_positions(phi):
    # Posisjonene til de tre pendelene i (x, y)-planet
    r1 = np.stack((L1*np.sin(phi[:,0]), -L1*np.cos(phi[:,0])), axis=1)
    r2 = r1 + np.stack((L2*np.sin(phi[:,1]), -L2*np.cos(phi[:,1])), axis=1)
    r3 = r2 + np.stack((L3*np.sin(phi[:,2]), -L3*np.cos(phi[:,2])), axis=1)
    return r1, r2, r3


# Pygame setup
pygame.init()
screen_width, screen_height = 800, 600
screen = pygame.display.set_mode((screen_width, screen_height))
clock = pygame.time.Clock()

# Colors
black = (0, 0, 0)
pendulum_colors = [(66, 124, 179), (66, 124, 179), (66, 124, 179)]

# Scale factor for drawing
scale_factor = 100

# Animation loop
running = True
current_frame = 0
while running:
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            running = False

    screen.fill(black)

    # Calculate positions for current frame
    R1, R2, R3 = r1[current_frame], r2[current_frame], r3[current_frame]

    pygame.draw.line(screen, pendulum_colors[0], (screen_width // 2, screen_height // 2), (screen_width // 2 + R1[0] * scale_factor, screen_height // 2 + R1[1] * scale_factor), 2)
    pygame.draw.line(screen, pendulum_colors[1], (screen_width // 2 + R1[0] * scale_factor, screen_height // 2 + R1[1] * scale_factor), (screen_width // 2 + R2[0] * scale_factor, screen_height // 2 + R2[1] * scale_factor), 2)
    pygame.draw.line(screen, pendulum_colors[2], (screen_width // 2 + R2[0] * scale_factor, screen_height // 2 + R2[1] * scale_factor), (screen_width // 2 + R3[0] * scale_factor, screen_height // 2 + R3[1] * scale_factor), 2)
    # Draw pendulums
    for i, (r, color) in enumerate(zip([R1, R2, R3], pendulum_colors)):
        pygame.draw.circle(screen, color, (int(screen_width // 2 + r[0] * scale_factor), int(screen_height // 2 + r[1] * scale_factor)), int(10 * ((3 * [m1, m2, m3][i]) / (4 * np.pi)) ** (1 / 3)))

    pygame.display.flip()
    current_frame += 1
    if current_frame >= N:
        current_frame = 0
    clock.tick(60)

pygame.quit()
sys.exit()
