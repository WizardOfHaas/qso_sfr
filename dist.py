import sys
import math
import scipy.integrate as integrate
import scipy.special as special
import numpy as np
import matplotlib.pyplot as plt

#Constant Definitions
#Note: most equations have already been numericalized with Mathematica
#http://planck.caltech.edu/pub/2015results/Planck_2015_Results_XIII_Cosmological_Parameters.pdf
Om = 0.308
OA = 1 - Om #Assume flat
Ok = 1 - Om - OA
c = 299792458
H0 = 1.62e-18
pc = 3.086e16
ly = 9.4607e15
L_sun = 3.828e26
M_sun = 4.83

#Equation definitions: distances
def E(z):
	return math.sqrt(OA + Om * (1 + z)**3)

def dH(z):
	return (c * z) / (H0 * E(z))

def dC(z):
	return dH(z) * integrate.quad(lambda z: 1 / E(z), 0, z)[0]

def dL(z):
	return (1 + z) * dC(z)

#Function definitions: m'bins
#Expects..
#	arr: 2d array; dependant, independant variable
#	dBin: bin size
#	fun: function to reduce to each bin

def bin_it(arr, dBin, fun):
	bins = []

	bMin = 0 #np.amin(arr[:,0])
	bMax = 8 #np.amax(arr[:,0])

	for b in np.arange(bMin, bMax, dBin):
		bin_arr = np.array(filter(lambda x: x[0] > b and x[0] <= b + dBin, arr))

		bin_mean = 0
		bin_std = 0
		bin_per = 0

		if len(bin_arr) > 0:
			bin_mean = np.mean(bin_arr[:,1])
			bin_std = np.std(bin_arr[:,1]) / bin_mean
			bin_per = len(bin_arr) / float(len(arr))

			#bin_arr = np.array(filter(lambda x: x[1] < 3 * bin_std + bin_mean, bin_arr))

			#if len(bin_arr) > 0:
			#	bin_val = fun(bin_arr[:,1])

		bins.append([
			b + dBin / 2,
			bin_mean,
			bin_std,
			bin_per
		])

	return np.array(bins)

#Function definition: wrapper to generate data
def calc_sfr(path):
	file = path
	csv_file = open(file, 'r')
	header = csv_file.readline().rstrip()

	sfr_data = []

	for l in csv_file:
		l = l.rstrip()
		data = l.split(",")

		z = float(data[6])
		f_U = float(data[5])

		if z > 0 and f_U > 0:
			d = dL(z)

			#From SDSS
			m_U = 22.5 - 2.5 * math.log(f_U, 10)

			#Derived from https://en.wikipedia.org/wiki/Luminosity
			L_U = (d / ly) ** 2 * 10 ** (-(m_U + 2.72) / 2.5) * L_sun

			#SFR_U = L_U * 10000000. * 3e-47 * 33545 #https://ned.ipac.caltech.edu/level5/Sept12/Calzetti/Calzetti1_2.html
			SFR_U = L_U * 5.9e-21 #http://www.physics.usyd.edu.au/~ahopkins/thesis/node11.html

			data.append(str(d))
			data.append(str(L_U))
			data.append(str(SFR_U))

			sfr_data.append([
				z,
				SFR_U
			])

	sfr_data = np.array(sfr_data)
	return sfr_data

qso_sfr_data = calc_sfr('data/qso_data.csv')
gal_sfr_data = calc_sfr('data/gal_data.csv')

dz = 0.5
bin_mean = lambda a: np.mean(a)
bin_std = lambda a: np.std(a)

qso_sfr_binned = bin_it(qso_sfr_data, dz, bin_mean)

gal_sfr_binned = bin_it(gal_sfr_data, dz, bin_mean)

plt.figure(1)

plt.subplot(311)

plt.errorbar(qso_sfr_binned[:,0], qso_sfr_binned[:,1], yerr=qso_sfr_binned[:,2], label='QSO', fmt="r-o")
#plt.plot(qso_sfr_binned[:,0], qso_sfr_binned[:,2], 'r-.')

plt.errorbar(gal_sfr_binned[:,0], gal_sfr_binned[:,1], yerr=gal_sfr_binned[:,2], label='GAL', fmt="b-o")
#plt.plot(gal_sfr_binned[:,0], gal_sfr_binned[:,2] , 'b-.')

plt.title("U Band Derived SFR Histories")
plt.legend(loc='upper left')
plt.xlabel("z")
plt.ylabel("SFR [log(M sun yr^-1)]")

plt.subplot(312)
plt.plot(qso_sfr_binned[:,0], (qso_sfr_binned[:,1] - gal_sfr_binned[:,1]) / gal_sfr_binned[:,1], 'g-o')

plt.title("SFR % Diff (QSO / GAL)")
plt.xlabel("z")
plt.ylabel("%")

plt.subplot(313)

plt.plot(qso_sfr_binned[:,0], qso_sfr_binned[:,3], 'r-o', label='QSO')
plt.plot(gal_sfr_binned[:,0], gal_sfr_binned[:,3], 'b-o', label='GAL')

plt.title("% Pop at Given Z")
plt.xlabel("z")
plt.ylabel("% of Pop")

plt.show()