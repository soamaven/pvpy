""" Test for photovoltaic.pv """
import numpy as np
import photovoltaic as pv
import matplotlib.pyplot as plt

Se = 100
Le = 1e-4
De = 4
xj = 1e-4

print('emitter SL/D',Se*Le/De)

xj_Wd = xj
Sr = 100
Wb = 0.01
Lb = 0.03
Db = 10
print('base SL/D',Sr*Lb/Db)
wavelength, abs_coeff, nd, kd = pv.optical_properties()
plt.plot(wavelength, abs_coeff)
plt.xlim(300,1200)
plt.ylim(1e-2,1e6)
plt.semilogy()
plt.show()

ab= abs_coeff

def IQE_baseL(ab, xj_Wd, Sr, Wb, Lb, Db):
    """return quantum efficiency of the base
    Given:
    ab - (cm) absorption coefficient
    xj_Wd - junction depth (cm)
    Sr - surface recombination velocity (cm/s)
    Lb - diffusion lenght of minority carrier in the base (cm)
    Db - diffusivity of minority carriers in the base (cm²/Vs)
    """
    GF = (ab * Lb - (
        (Sr * Lb / Db) * (np.cosh(Wb / Lb) - np.exp(-ab * Wb)) + np.sinh(Wb / Lb) + Lb * ab * np.exp(-ab * Wb)) / (
              (Sr * Lb / Db) * np.sinh(Wb / Lb) + np.cosh(Wb / Lb)))
    QEB = (np.exp(-ab * (xj_Wd))) * (Lb * ab / (ab**2 * Lb**2 - 1)) * GF
    return QEB, GF

IQEE = pv.IQE_emitter(abs_coeff, Se, Le, De, xj)
IQEB, GF = IQE_baseL(abs_coeff, xj_Wd, Sr, Wb, Lb, Db)
(Lb * ab / (ab**2 * Lb**2 - 1))
IQEt = IQEE+IQEB

plt.plot(wavelength, (np.exp(-ab * (xj_Wd))) * (Lb * ab / (ab**2 * Lb**2 - 1)))
plt.semilogy()

plt.show()

plt.plot(wavelength, IQEE, label='emitter')
plt.plot(wavelength, IQEB, label='base')
plt.plot(wavelength, IQEt, label='IQE')
plt.xlim(300,1200)
plt.ylim(0,1)
plt.legend(loc=0)
plt.show()
