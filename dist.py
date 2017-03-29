import sys
import math
import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import matplotlib.pyplot as plt

#Constant Definitions
#Note: most equations have already been numericalized with Mathematica
Om = 0.27
OA = 0.73
Ok = 1 - Om - OA
c = 299792458
H0 = 69.3
pc = 3.086e16
ly = 1.057e-16
L_sun = 3.828e26
M_sun = 4.83

#Equation definitions: distances
def E(z):
	return math.sqrt(0.73 + 0.27 * (1 + z)**3)

def dH(z):
	return (c * z) / H0

def dC(z):
	return dH(z) * integrate.quad(lambda z: 1 / E(z), 0, z)[0]

def dL(z):
	return (1 + z) * dC(z)

#Function definitions: m'bins
#Expects..
#	arr: 2d array; dependant, independant variable
#	dBin: bin size
def bin_it(arr, dBin):
	bins = []

	bMin = np.amin(arr[:,0])
	bMax = np.amax(arr[:,0])

	for b in np.arange(bMin, bMax, dBin):
		bin_arr = np.array(filter(lambda x: x[1] > b and x[1] <= b + dBin, arr))
		bins.append([
			b,
			np.mean(bin_arr[:,0])
		])

	return np.array(bins)

file = sys.argv[1]
csv_file = open(file, 'r')
header = csv_file.readline().rstrip()
#header = header + ",d,L_U,SFR_U"
#print header

sfr_data = []

for l in csv_file:
	l = l.rstrip()
	data = l.split(",")

	z = float(data[6])
	f_U = float(data[5])

	if z > 0 and f_U > 0:
		d = dL(z)

		m_U = 22.5 - 2.5 * math.log(f_U, 10)
		L_U = (d / ly) ** 2 * 10 ** (-(m_U + 2.72) / 2.5) * L_sun

		SFR_U = L_U * 3e-47 * 3.543e-7 #https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html

		data.append(str(d))
		data.append(str(L_U))
		data.append(str(SFR_U))

		sfr_data.append([
			z,
			SFR_U
		])

		#print ",".join(data)

sfr_data = np.array(sfr_data)

sfr_binned = bin_it(sfr_data, 0.5)

plt.plot(sfr_binned[:,0], sfr_binned[:,1], 'ro', label='Z Binned SFR')
plt.show()