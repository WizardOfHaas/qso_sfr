import math
import matplotlib.pyplot as plt
import numpy as np

qso_flux_file = open('qso_flux.csv', 'r')
qso_flux = []

L_sun = 3.828e26
M_sun = 4.83
c = 299792458
H0 = 2.2e-18
pc = 3.086e16

#Read in and cconvert to SFR
for l in qso_flux_file:
	data = l.split(",")
	
	z = float(data[5].rstrip())
	f_U = float(data[3]) * (10 ** -9)

	if f_U > 0 and z > 0:
		d = ((z * c) / H0) / pc

		m_U = 22.5 - 2.5 * math.log(f_U, 10)
		M_U = 2.5 * math.log(d, 10) + 5 + m_U
		L_U = 10 ** ((M_sun - M_U)/2.5) * L_sun

		qso_flux.append([
			L_U / 5.9e21,
			z
		])

#Bin data by z
sfr_mean = []
dz = 0.5
for z in np.arange(0, 7, dz):
	sfr_mean.append(
		[
			np.mean(filter(lambda x: x[1] > z and x[1] <= z + dz, qso_flux)),
			z
		]
	)

sfr_mean = np.array(sfr_mean)

plt.plot(sfr_mean[:,1], sfr_mean[:,0], 'o')
plt.show()