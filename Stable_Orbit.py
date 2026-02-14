"""
Newtonian Two-Body Gravity Simulation (RK4)
-------------------------------------------
Author: [Bhaskar Chatterjee]
Date: Sept 2025
Description: 
    A numerical simulation of the Keplerian two-body problem using a 
    4th-Order Runge-Kutta (RK4) integrator. This script demonstrates 
    orbital stability and energy conservation in a Newtonian potential.
"""
import numpy as np
import matplotlib.pyplot as plt 

# Define toy values for G, m1, m2
G= 1.0
m_1=2.0
m_2=1.0
"""
    4th-Order Runge-Kutta Integrator (RK4).
    
    Parameters:
        fun : callable, Right-hand side of the ODE (dy/dt = fun(t, y))
        dt  : float, Time step
        y0  : array, Initial state vector
        t0  : float, Start time
        tf  : float, End time
        
    Returns:
        t : array, Time grid
        y : array, Solution state vector at each time step
    """
def rk_4(fun,dt,y0,t0,tf):
    t=np.arange(t0,tf+dt,dt)
    n=len(t)
    y=np.zeros((n,len(y0)))
    y[0]=y0
    for i in range(n-1):
        k1=fun(t[i],y[i])
        k2=fun(t[i]+dt/2,y[i]+k1*dt/2)
        k3=fun(t[i]+dt/2,y[i]+k2*dt/2)
        k4=fun(t[i]+dt,y[i]+k3*dt)
        y[i+1]=y[i]+(k1+2*k2+2*k3+k4)*dt/6
    return t,y
# Mathematical function for the two body problem
def two_body(t,y):
    x1,y1,x2,y2,vx1,vy1,vx2,vy2=y
    # Calculate separation distance r
    # Note: Avoid singularity at r=0 in real-world code
    r=np.sqrt((x2-x1)**2+(y2-y1)**2)
    ax1=G*m_2*(x2-x1)/r**3
    ay1=G*m_2*(y2-y1)/r**3
    ax2=G*m_1*(x1-x2)/r**3
    ay2=G*m_1*(y1-y2)/r**3
    return np.array([vx1,vy1,vx2,vy2,ax1,ay1,ax2,ay2])
# Initial conditions
dt=0.0001     # Small time step for high accuracy
t0=0.0
tf=20.0       # Duration of the simulation
y0=np.array([1.0,0.0,-1.0,0.0,0.0,0.25,0.0,-0.5]) #in the form [x1,y1,x2,y2,vx1,vy1,vx2,vy2]

# Run the RK4 method to solve the two body problem
t,y=rk_4(two_body,dt,y0,t0,tf)
# Extract positions for plotting as they are stored in y as a 2D array with columns corresponding to x1,y1,x2,y2,vx1,vy1,vx2,vy2
x1=y[:,0]
y1=y[:,1]
x2=y[:,2]
y2=y[:,3]   
# Calculate the total energy of the system at each time step to check for conservation of energy (should be constant for a closed system)
def energy(t,y):
    x1,y1,x2,y2,vx1,vy1,vx2,vy2=y.T
    r=np.sqrt((x2-x1)**2+(y2-y1)**2)
    K=0.5*m_1*(vx1**2+vy1**2)+0.5*m_2*(vx2**2+vy2**2)
    U=-G*m_1*m_2/r
    return K+U
E=energy(t,y)
#Plotting the trajectories of the two bodies and the total energy of the system
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
#Orbit
ax1.plot(x1,y1,label='Body 1')
ax1.plot(x2,y2,label='Body 2')
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.set_title('Two Body Problem using RK4')
ax1.legend()
ax1.axis('equal')
ax1.grid()

 #Total Energy
ax2.plot(t,E)
ax2.set_xlabel('Time') 
ax2.set_ylabel('Total Energy')
ax2.set_title('Total Energy of the Two Body System')
ax2.grid()
plt.show()         