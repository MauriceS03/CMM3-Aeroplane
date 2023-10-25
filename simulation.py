'''

# importing modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
import constants
import Forms

#-----------------------------------------------------------------------------------------------------------
# User parameters

velocity = 100  # Aircraft velocity in m/s
gamma = 0.00   # Path angle in rad

gravity = 9.81  # Gravitational acceleration in m/s^2
air_density = 1.0065    # Air density in kg/m^3
wing_surface = 20.0 # Wing surface in m^2
cbar = 1.75 # airfoil chord in m
mass = 1300.0     # Mass of the airplane in kg
inertia_yy = 7000   # Moment of inertia in kg/m^2

#-----------------------------------------------------------------------------------------------------------
# Defining functions for later use. The functions are derived from the project breif.

def CL (a, d):
    return Forms.Coefficient_of_Lift(a, d)

def CM (a, d):
    return Forms.Coefficient_of_Moment(a, d)

def CD (a, d):
    return Forms.Coefficient_of_Drag(a, d)

def L (a, d):
    return Forms.Lift(a, d)

def D (a, d):
    return Forms.Drag(a, d)

def M (a, d):
    return Forms.Moment(a, d)

def Thrust (a, d, the):
    return Forms.Engine_Thrust(a, d, the)

# Runge-Kutta method for integral solving. Parameters from dx/dt = f
def RK4(x, f, dt):
    k1 = f
    k2 = f + 0.5 * dt * k1
    k3 = f + 0.5 * dt * k2
    k4 = f + dt * k3
    
    return x + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

#-----------------------------------------------------------------------------------------------------------
# Final calculations & output

# Define the equilibrium equation as f(a)
def f(a):
    return Forms.Equilibrium(a)   

# Solve for alpha and delta
initial_guess = 0.01  # Provide an initial guess
alpha = newton(f, initial_guess)
delta = -(constants.CM0 + constants.CMa * alpha)/constants.CMde

# Calculating other variables to output
theta = alpha + gamma
ub = velocity * np.cos(alpha)
wb = velocity * np.sin(alpha)

# Calculating Thrust
thrust = Thrust(alpha, delta, theta) 

print(f"alpha = {alpha}")
print(f"delta = {delta}")
print(f"thrust = {thrust}")
print(f"theta = {theta}")
print(f"ub = {ub}")
print(f"wb = {wb}")
#Moment replaced with M
print(M(alpha, delta))

#-----------------------------------------------------------------------------------------------------------
# Applying Euler method for solving differential DOF equations

t0 = 0      # Initial time (s)
tEnd = 10   # End time (s)
dt = 0.1    # Time step size (s)

q = 0 # Initial angular velocity in rad/sec
xe = 0
ze = 0
t = t0

moment = 0

tValues = [t0]
thetaValues = [theta]
qValues = [q]
xeValues = [xe]
zeValues = [ze]
ubValues = [ub]
wbValues = [wb]
gammaValues = [gamma]
alphaValues = [alpha]
momentValues = [moment]

while t < tEnd:
    if t >= 1:
        delta = -0.0572
    # Compute new values using the DOF equations
    theta = RK4(theta, q, dt)
    alpha = np.arctan2(wb, ub)
    gamma = theta - alpha
    #Moment replaced with M
    moment = M(alpha, delta)
    thrust = Thrust(alpha, delta, theta)
    q = RK4(q, (moment/inertia_yy), dt)
    xe += (ub * np.cos(theta) + wb * np.sin(theta)) * dt
    ze -= (- ub * np.sin(theta) + wb * np.cos(theta)) * dt
    ub += (L(alpha, delta) * np.sin(alpha) / mass - D(alpha, delta) *
           np.cos(alpha) / mass - q * wb - gravity * np.sin(theta) +
           thrust/mass) * dt
    wb += (- L(alpha, delta) * np.cos(alpha) / mass - D(alpha, delta) *
           np.sin(alpha) / mass + q * ub + gravity * np.cos(theta)) * dt
    # Append new values to arrays
    t += dt
    tValues.append(round(t, 1))
    thetaValues.append(theta)
    qValues.append(q)
    xeValues.append(xe)
    zeValues.append(ze)
    ubValues.append(ub)
    wbValues.append(wb)
    alphaValues.append(alpha)
    gammaValues.append(gamma)
    momentValues.append(moment)
    '''
    plt.subplot(4, 2, 1)
    plt.plot(tValues, ubValues)
    plt.xlabel('time')
    plt.ylabel('ub')
    plt.subplot(4, 2, 2)
    plt.plot(tValues, wbValues)
    plt.xlabel('time')
    plt.ylabel('wb')
    plt.subplot(4, 2, 3)
    plt.plot(tValues, qValues)
    plt.xlabel('time')
    plt.ylabel('q')
    plt.subplot(4, 2, 4)
    plt.plot(tValues, thetaValues)
    plt.xlabel('time')
    plt.ylabel('theta')
    plt.subplot(4, 2, 5)
    plt.plot(tValues, gammaValues)
    plt.xlabel('time')
    plt.ylabel('path angle')
    plt.subplot(4, 2, 6)
    plt.plot(tValues, zeValues)
    plt.xlabel('time')
    plt.ylabel('ze')
    plt.subplot(4, 2, 7)
    plt.plot(tValues, alphaValues)
    plt.xlabel('time')
    plt.ylabel('alpha')
    plt.subplot(4, 2, 8)
    plt.plot(tValues, momentValues)
    plt.xlabel('time')
    plt.ylabel('moment')
    #plt.tight_layout()
    plt.show()
    
    plt.subplot(2, 2, 1)
    plt.plot(tValues, momentValues)
    plt.ylabel('moment')
    plt.grid()
    plt.subplot(2, 2, 2)
    plt.plot(tValues, qValues)
    plt.ylabel('q')
    plt.grid()
    plt.subplot(2, 2, 3)
    plt.plot(tValues, alphaValues)
    plt.ylabel('alpha')
    plt.grid()
    plt.subplot(2, 2, 4)
    plt.plot(tValues, thetaValues)
    plt.ylabel('theta')
    plt.grid()
    plt.tight_layout()
    plt.show()
    print(moment, q)
    print(tValues)
    '''
    #print(delta)

print(f"new alpha = {alpha}")
print(f"new delta = {delta}")
print(f"new ub = {ub}")
print(f"new wb = {wb}")
print(f"new q = {q}")
print(f"new gamma = {gamma}")
print(f"new theta = {theta}")


# Plot the results
#plt.plot(tValues, alphaValues, 'b-')
plt.subplot(4, 2, 1)
plt.plot(tValues, ubValues)
plt.xlabel('time')
plt.ylabel('ub')
plt.subplot(4, 2, 2)
plt.plot(tValues, wbValues)
plt.xlabel('time')
plt.ylabel('wb')
plt.subplot(4, 2, 3)
plt.plot(tValues, qValues)
plt.xlabel('time')
plt.ylabel('q')
plt.subplot(4, 2, 4)
plt.plot(tValues, thetaValues)
plt.xlabel('time')
plt.ylabel('theta')
plt.subplot(4, 2, 5)
plt.plot(tValues, gammaValues)
plt.xlabel('time')
plt.ylabel('path angle')
plt.subplot(4, 2, 6)
plt.plot(tValues, zeValues)
plt.xlabel('time')
plt.ylabel('ze')
plt.subplot(4, 2, 7)
plt.plot(tValues, alphaValues)
plt.xlabel('time')
plt.ylabel('alpha')
plt.subplot(4, 2, 8)
plt.plot(tValues, momentValues)
plt.xlabel('time')
plt.ylabel('moment')
#plt.tight_layout()
plt.show()
