import numpy as np
import math
import matplotlib.pyplot as plt

c = 299792458 # speed of line in m/s
G = 6.6743E-11 # Newton's gravitational constant in m^3/kg*s^2
rho_c = 9.9E17 # central density in kg/m^3

R_NS = 11400 # Radius of neutron star in m
M_NS = 1.4*1.9891E30 # Mass of neutron star in kg
rho_c_tol = (15.0*M_NS)/(8.0*math.pi*R_NS**3)

r = np.linspace(0, R_NS, 1000)
xi = r/R_NS
C = ((8.0*math.pi/15.0)*R_NS**2*rho_c_tol)*(G/c**2)
print(C)

rho_tol = (15.0*C*(1-xi**2))/(8.0*math.pi*R_NS**2)
m_tol = R_NS*C*((5.0/2.0)*xi**3 - (3.0/2.0)*xi**5)
C_1 = 1.0 - 5.0*C/3.0
C_2 = math.atan(math.sqrt(C/(3.0*(1.0 - 2.0*C)))) + 0.5*math.log((1.0/6.0) + math.sqrt((1.0-2.0*C)/(3.0*C)))
E_lamb_tol = 1.0 - C*xi**2*(5.0 - 3.0*xi**2)
phi_tol = C_2 - 0.5*np.log(xi**2 - (5.0/6.0) + np.sqrt(E_lamb_tol/(3.0*C)))
E_nu_tol = C_1*np.cos(phi_tol)*np.cos(phi_tol)

p_tol = ((np.sqrt(3.0*C*E_lamb_tol)*np.tan(phi_tol) - C*(5.0 - 3.0*xi**2)/2.0)/(4.0*math.pi*R_NS**2))

p_c = 0.9*p_tol[0]

plt.plot(r/R_NS, p_tol/p_c)
plt.show()
