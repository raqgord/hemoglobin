# This script fits oxygen hemoglobin binding data with the Adair function and the Pauling function.

# Written by Raquel Gordon (2017). Last modified January 2017.
# hemoglobinFerryAndGreenExp.py is licensed under a 
# Creative Commons Attribution-ShareAlike 3.0 Unported License.


import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
from scipy.optimize import curve_fit

def funcAdair(p, K):
	K1 = K
	K2 = K*2.8
	K3 = K*14.0
	K4 = K*160.0
	# First pass affinity constant ratios determined by trail and error.
	# Adair equation as described in L. Pauling, PNAS (1935)
	return (((K1 * p) + (2.0 * K2 * p**2) + (3.0 * K3 * p**3) + (4.0 * K4 * p**4)) / \
			(4.0*(1.0 + (K1 * p) + (K2 * p**2) + (K3 * p**3) + (K4 * p**4))))

xdatapH738 = np.array([3.7, 5.5, 7.7, 8.1, 9.9, 11.2, 12.2, 14.8, 18.8, 27.0, 28.6])
print np.size(xdatapH738)
xdata = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]

ydatapH738 = np.array([0.126, 0.219, 0.433, 0.423, 0.558, 0.7, 0.697, 0.805, 0.888, 0.949, 0.956])
print np.size(ydatapH738)
popt, pcov = curve_fit(funcAdair, xdatapH738, ydatapH738)
print popt
yFit = funcAdair(xdatapH738, popt)

a1 = plt.subplot(211)
a1.plot(xdatapH738, ydatapH738, 'r-', label='Ferry and Green Data')
a1.plot(xdatapH738, yFit, 'g-', label='Adair Fit')
legend = a1.legend(loc='lower right')

resid = yFit - ydatapH738
plt.subplot(212)
plt.plot(xdatapH738, resid, 'y-')
plt.ylabel('Residual Fraction Bound')
plt.xlabel('Residual Oxygen')

plt.ylabel('Fraction Bound')
plt.xlabel('Oxygen Pressure')
plt.title('Hemoglobin Adair Fit')
plt.show()


#Fit with Pauling equation
def funcPauling(p, K, alpha):
	return ((K * p) + ((2 * alpha + 1) * K**2 * p**2) + (3 * alpha**2 * K**3 * p**3) + (alpha**4 * K**4 * p**4)) / \
			(1 + (4 * K * p) + ((4*alpha + 2) * K**2 * p**2) + (4 * alpha**2 * K**3 * p**3) + (alpha**4 * K**4 * p**4))

popt, pcov = curve_fit(funcPauling, xdatapH738, ydatapH738)
print popt
yFit = funcPauling(xdatapH738, popt[0], popt[1])

xdatapH738 = np.array([3.7, 5.5, 7.7, 8.1, 9.9, 11.2, 12.2, 14.8, 18.8, 27.0, 28.6])
ydatapH738 = np.array([0.126, 0.219, 0.433, 0.423, 0.558, 0.7, 0.697, 0.805, 0.888, 0.949, 0.956])

resid = yFit - ydatapH738
plt.subplot(212)
plt.plot(xdatapH738, resid, 'y-')
plt.ylabel('Residual Fraction Bound')
plt.xlabel('Residual Oxygen')

a2 = plt.subplot(211)
a2.plot(xdatapH738, ydatapH738, 'r-', label='Ferry and Green Data')
a2.plot(xdatapH738, yFit, 'b-', label='Pauling Fit')
legend = a2.legend(loc='lower right')
plt.ylabel('Fraction Bound')
plt.xlabel('Oxygen Pressure')
plt.title('Hemoglobin Pauling Fit')
plt.show()
