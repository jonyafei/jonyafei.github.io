import numpy as np
import sys


kb = 1.38065E-23  # Boltzmann constant, unit: J/K
h = 6.626069E-34  # Planck constant, unit: J*s
R = 8.3145  # Gas constant J/(mol*K)
kcal2J = 4.184043E+03  # 1 kcal/mol = 4.184043 kJ/mol
eV2J = 96.48534E+03  # 1eV = 96.48534 kJ/mol

deltaGa = float(sys.argv[1])  # Activation free energy, unit: kcal/mol
T = float(sys.argv[2])  # Temperature: K

rate = np.exp(- (deltaGa * kcal2J) / (R * T)) * (kb * T / h)
t_half = np.log(2) / rate

print('Your input activation free energy barrier is %.2f kcal/mol' % deltaGa)
print("The rate constant under %.2f K is %.3e s^-1" % (T, rate))
print("The half life-time of the first order reaction is %.3e s" % t_half)

# usage: python G2k.py [deltaGa] [Temperature]
# such as: python G2k.py 23.0 298.15
# unit of deltaG is kcal/mol
# unit of Temperature is K

