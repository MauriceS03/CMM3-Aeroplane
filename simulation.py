
# importing modules
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
import constants
import forms

#-----------------------------------------------------------------------------------------------------------
# User parameters

velocity = 100  # Aircraft velocity in m/s
gamma = 0.00   # Path angle in rad

#-----------------------------------------------------------------------------------------------------------
# Defining functions for later use. The functions are derived from the project breif.

def CL (a, d):
    return forms.Coefficient_of_Lift(a, d)

def CM (a, d):
    return forms.Coefficient_of_Moment(a, d)

def CD (a, d):
    return forms.Coefficient_of_Drag(a, d)

def L (a, d):
    return forms.Lift(a, d, velocity, gamma)

def D (a, d):
    return forms.Drag(a, d, velocity, gamma)

def M (a, d):
    return forms.Moment(a, d, velocity, gamma)

def Thrust (a, d, the):
    return forms.Engine_Thrust(a, d, the, velocity, gamma)

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
    return forms.Equilibrium(a, velocity, gamma)   

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
    q = RK4(q, (moment/constants.inertia_yy), dt)
    xe += (ub * np.cos(theta) + wb * np.sin(theta)) * dt
    ze -= (- ub * np.sin(theta) + wb * np.cos(theta)) * dt
    ub += (L(alpha, delta) * np.sin(alpha) / constants.mass - D(alpha, delta) *
           np.cos(alpha) / constants.mass - q * wb - constants.gravity * np.sin(theta) +
           thrust/constants.mass) * dt
    wb += (- L(alpha, delta) * np.cos(alpha) / constants.mass - D(alpha, delta) *
           np.sin(alpha) / constants.mass + q * ub + constants.gravity * np.cos(theta)) * dt
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

# Plot the results
plt.plot(tValues, alphaValues, 'b-')
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
plt.show()
