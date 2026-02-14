"""
Binary Black Hole Inspiral Simulation (Dissipative Newtonian)
-------------------------------------------------------------
Author: [Bhakar Chatterjee]
Date: Oct 2026
Context: Gravitational Wave Physics Preparation

Description:
    This script simulates the 'Inspiral' phase of a binary black hole merger.
    It extends the standard Newtonian Two-Body problem by adding a 
    phenomenological dissipative term (radiation reaction) to the equations of motion.

    The simulation demonstrates:
    1. Orbital Decay: The transition from stable orbit to merger.
    2. Energy Loss: The system loses mechanical energy over time (mimicking GW emission).
    
    System Parameters (Toy Units):
    - G = 1.0 (Gravitational Constant)
    - m1 = 10.0, m2 = 10.0 (Equal Mass Binary)
    - gamma = 0.05 (Dissipation Strength)
"""
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt 
from matplotlib.animation import FuncAnimation

# --- PHYSICAL CONSTANTS ---
G = 1.0   
m_1 = 10.0  # Heavier masses for stronger gravity
m_2 = 10.0
gamma = 0.04 # DISSIPATION FACTOR: Controls the rate of inspiral (Energy Loss)

def rk_4(fun, dt, y0, t0, tf):
    """
    4th-Order Runge-Kutta Integrator (Standard Solver).
    """
    t = np.arange(t0, tf + dt, dt)
    n = len(t)
    y = np.zeros((n, len(y0)))
    y[0] = y0
    
    for i in range(n - 1):
        k1 = fun(t[i], y[i])
        k2 = fun(t[i] + dt/2, y[i] + k1 * dt/2)
        k3 = fun(t[i] + dt/2, y[i] + k2 * dt/2)
        k4 = fun(t[i] + dt, y[i] + k3 * dt)
        y[i+1] = y[i] + (k1 + 2*k2 + 2*k3 + k4) * dt / 6
        
    return t, y

def inspiral_dynamics(t, y):
    """
    Equations of Motion with Dissipation.
    a_total = a_gravity + a_dissipative
    """
    x1, y1, x2, y2, vx1, vy1, vx2, vy2 = y
    
    rx = x2 - x1
    ry = y2 - y1
    r = np.sqrt(rx**2 + ry**2)
    
    if r < 0.1: # Merger condition (stop singularity)
        return np.zeros(8) 

    # 1. NEWTONIAN GRAVITY (The Pull) 
    # F = Gm1m2/r^2
    ax1_grav = G * m_2 * (rx) / r**3
    ay1_grav = G * m_2 * (ry) / r**3
    ax2_grav = G * m_1 * (-rx) / r**3
    ay2_grav = G * m_1 * (-ry) / r**3
    
    # 2. DISSIPATIVE TERM (The Drag) 
    # We add a term proportional to velocity to mimic radiation reaction.
    # In full GR, this depends on the quadrupole moment derivative.
    # Here, we model it as: a_diss = -gamma * v * (1/r^2)
    # The (1/r^2) factor makes the inspiral speed up as they get closer (Chirp).
    
    scaling = 1.0 / (r**2) # Stronger drag at close range
    
    ax1_diss = -gamma * vx1 * scaling
    ay1_diss = -gamma * vy1 * scaling
    ax2_diss = -gamma * vx2 * scaling
    ay2_diss = -gamma * vy2 * scaling
    
    # Total Acceleration
    ax1 = ax1_grav + ax1_diss
    ay1 = ay1_grav + ay1_diss
    ax2 = ax2_grav + ax2_diss
    ay2 = ay2_grav + ay2_diss
    
    return np.array([vx1, vy1, vx2, vy2, ax1, ay1, ax2, ay2])

# --- INITIAL CONDITIONS ---
dt = 0.002       
t0 = 0.0
tf = 60.0        # Longer time to see the spiral

# High Eccentricity Setup:
# Start far apart (d=4.0) with low tangential velocity (v=0.8).
# This creates a "falling" orbit rather than a circular one.
y0 = np.array([2.0, 0.0, -2.0, 0.0, 0.0, 0.8, 0.0, -0.8]) 

# --- RUN SIMULATION ---
print(f"Simulating Inspiral (gamma={gamma})...")
t, y = rk_4(inspiral_dynamics, dt, y0, t0, tf)

# Extract trajectories
x1, y1 = y[:, 0], y[:, 1]
x2, y2 = y[:, 2], y[:, 3]
vx1, vy1 = y[:, 4], y[:, 5]
vx2, vy2 = y[:, 6], y[:, 7]

# --- ENERGY ANALYSIS ---
def get_energy(t, y_state):
    """Calculates Total Mechanical Energy (K + U)."""
    x1, y1, x2, y2, vx1, vy1, vx2, vy2 = y_state.T
    r = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
    
    # Filter out merged steps (r=0) to avoid plotting errors
    mask = r > 0.5
    
    K = 0.5 * m_1 * (vx1**2 + vy1**2) + 0.5 * m_2 * (vx2**2 + vy2**2)
    U = -G * m_1 * m_2 / r
    return K + U, mask

E, mask = get_energy(t, y)

# --- VISUALIZATION ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# --- PLOT 1: The Orbit ---
ax1.plot(x1, y1, color='purple', alpha=0.3, linewidth=1)
ax1.plot(x1, y1, label='Black Hole 1', color='indigo', linewidth=1.5)
ax1.plot(x2, y2, label='Black Hole 2', color='darkorange', linewidth=1.5)
ax1.set_xlabel('x ($R_s$)')
ax1.set_ylabel('y ($R_s$)')
ax1.set_title('Eccentric Binary Inspiral')
ax1.legend(loc='upper right')
ax1.axis('equal')
ax1.grid(True, alpha=0.2)

# --- PLOT 2: The Energy (Full View) ---
r_array = np.sqrt((x2 - x1)**2 + (y2 - y1)**2)
# Filter out the crash for clean plotting (removes singularity noise)
clean_mask = r_array > 0.1
t_clean = t[clean_mask]
E_clean = E[clean_mask]

ax2.plot(t_clean, E_clean, color='green', linewidth=2)
ax2.set_xlabel('Time (s)') 
ax2.set_ylabel('Total Energy')
ax2.set_title('Gravitational Wave Emission')
ax2.grid(True, alpha=0.2)
ax2.annotate('Merger Plunge', 
             xy=(t_clean[-1], E_clean[-1]), 
             xytext=(t_clean[-1]-15, E_clean[-1]+50),
             arrowprops=dict(facecolor='black', shrink=0.05),
             fontsize=10, fontweight='bold')

#  (The Zoom on the Stairs) 
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

# Create a box inside the main plot (width=40%, height=40%)
ax_ins = inset_axes(ax2, width="40%", height="40%", loc='lower left', borderpad=3)

# Plot the same data again, but zoomed in
ax_ins.plot(t_clean, E_clean, color='darkgreen', linewidth=1.5)

# Set the ZOOM limits (0 to 35s, and just the top of the energy)
ax_ins.set_xlim(0, 20)
ax_ins.set_ylim(E_clean[0]-10, E_clean[0]+0.5) # Tuning y-limits to see steps

ax_ins.set_title("Zoom: Pericenter Steps", fontsize=9)
ax_ins.grid(True, alpha=0.3)
ax_ins.tick_params(labelsize=8)

# Highlight the area we are zooming in on (Optional, adds lines connecting box)
from matplotlib.patches import Rectangle
rect = Rectangle((0, E_clean[0]-10), 35, 5.5, linewidth=1, edgecolor='gray', facecolor='none', linestyle='--')
ax2.add_patch(rect)

plt.tight_layout()
plt.savefig('binary_inspiral_inset.png', dpi=300)
print("Saved plot with Inset Zoom!")

# Setup the figure for animation
fig_anim, ax_anim = plt.subplots(figsize=(8, 8))
ax_anim.set_xlim(-2.5, 2.5)
ax_anim.set_ylim(-2.5, 2.5)
ax_anim.set_title(f"Binary Black Hole Inspiral (gamma={gamma})")
ax_anim.set_xlabel("x (m)")
ax_anim.set_ylabel("y (m)")
ax_anim.grid(True, alpha=0.3)
ax_anim.set_aspect('equal')

# Create the empty line objects (Bodies and Trails)
body1, = ax_anim.plot([], [], 'o', color='purple', markersize=12, label='BH 1')
body2, = ax_anim.plot([], [], 'o', color='orange', markersize=12, label='BH 2')
trail1, = ax_anim.plot([], [], '-', color='purple', linewidth=1, alpha=0.5)
trail2, = ax_anim.plot([], [], '-', color='orange', linewidth=1, alpha=0.5)
ax_anim.legend(loc='upper right')

# Initialization function: Plot the background of each frame
def init():
    body1.set_data([], [])
    body2.set_data([], [])
    trail1.set_data([], [])
    trail2.set_data([], [])
    return body1, body2, trail1, trail2

# Animation function: Update the plot for each frame
# We skip frames (step=10) to make the animation faster
skip = 20
def update(frame):
    i = frame * skip
    if i >= len(x1): return body1, body2, trail1, trail2
    
    # Update Body positions
    body1.set_data([x1[i]], [y1[i]])
    body2.set_data([x2[i]], [y2[i]])
    
    # Update Trails (History)
    # Show the entire history up to the current frame
    trail1.set_data(x1[:i], y1[:i])
    trail2.set_data(x2[:i], y2[:i])
    
    return body1, body2, trail1, trail2

# Create the animation object
frames = len(t) // skip
ani = FuncAnimation(fig_anim, update, frames=frames, init_func=init, blit=True, interval=20)

# SAVE THE GIF
ani.save('inspiral_movie.gif', writer='ffmpeg', fps=30)
print("Animation saved as 'inspiral_movie.gif'")

