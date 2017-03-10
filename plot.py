import math
import matplotlib.pyplot as plt
import numpy as np

#qso_flux_file = open('qso_flux.csv', 'r')
#qso_flux = []

L_sun = 3.828e26
M_sun = 4.83
c = 299792458
H0 = 2.2e-18
pc = 3.086e16

def sfr(file):
	#Read in and cconvert to SFR
	d_file = open(file, 'r')
	flux = []	

	for l in d_file:
		data = l.split(",")
	
		z = float(data[5].rstrip())
		f_U = float(data[3]) * (10 ** -9)

		if f_U > 0 and z > 0:
			d = ((z * c) / H0) / pc

			m_U = 22.5 - 2.5 * math.log(f_U, 10)
			M_U = 2.5 * math.log(d, 10) + 5 + m_U
			L_U = 10 ** ((M_sun - M_U)/2.5) * L_sun

			flux.append([
				L_U / 5.9e21,
				z
			])

	num_samples = len(flux) * 1.0

	#Bin data by z
	sfr_mean = []
	dz = 0.5
	for z in np.arange(0, 7, dz):
		z_bin = np.array(filter(lambda x: x[1] > z and x[1] <= z + dz, flux))

		sfr_mean.append(
			[
				math.log(np.mean(z_bin[:,0]), 10) if len(z_bin) > 0 else 0,
				z,
				len(z_bin) / num_samples
			]
		)
	sfr_mean = np.array(sfr_mean)

	#Bin data by SFR
	sfr_bin = []
	dsfr = 0.5
	for s in np.arange(0, 3.5, dsfr):
		s_bin = np.array(filter(lambda x: x[0] > s and x[0] <= s + dsfr, flux))

		sfr_bin.append(
			[
				s,
				len(s_bin) / num_samples
			]
		)

	sfr_bin = np.array(sfr_bin)
	
	return sfr_mean, sfr_bin

qso_sfr_z, qso_sfr_bin = sfr('qso_flux.csv')
gal_sfr_z, gal_sfr_bin = sfr('gal_flux.csv')

#print qso_sfr_z
#print qso_sfr_bin

plt.figure(1)

plt.subplot(311)
plt.plot(qso_sfr_z[:,1], qso_sfr_z[:,0], 'ro', label='QSO')
plt.plot(gal_sfr_z[:,1], gal_sfr_z[:,0], 'bo', label='Galaxy')

plt.title("U Band Derived SFR Histories: SDSS DR12")
plt.legend(loc='upper left')
plt.xlabel("z")
plt.ylabel("SFR (log(M sun yr^-1))")

plt.subplot(312)
plt.plot(qso_sfr_z[:,1], qso_sfr_z[:,2], 'ro', label='QSO')
plt.plot(gal_sfr_z[:,1], gal_sfr_z[:,2], 'bo', label='Galaxy')

plt.title("% Population at Given Z: SDSS DR12")
plt.xlabel("z")
plt.ylabel("Pop (%)")

plt.subplot(313)
plt.plot(qso_sfr_bin[:,0], qso_sfr_bin[:,1], 'ro', label='QSO')
plt.plot(gal_sfr_bin[:,0], gal_sfr_bin[:,1], 'bo', label='Galaxy')

plt.title("% Population with Given SFR: SDSS DR12")
plt.xlabel("SFR (M sun yr^-1)")
plt.ylabel("Pop (%)")

plt.show()
