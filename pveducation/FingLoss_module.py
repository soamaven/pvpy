import matplotlib.pyplot as plt
import numpy as np

import pveducation as pv

losses = pv.Pf_total(2, 0.035, 0.2, 1.6e-6, 50, 0.01, 0.001, 0.650)
Sf = np.linspace(0.01,0.5)
Ptotal, Presistivity, Pshading, Psheet  = pv.Pf_total(2.0,0.035,Sf, 1.6e-6, 50, 100e-4, 0.5e-4, 0.650)
plt.ylim(0,20)
plt.xlabel('finger spacing (cm)')
plt.ylabel('power loss (%)')
plt.plot(Sf, Ptotal, label='total')
plt.plot(Sf, Presistivity, label='resistive')
plt.plot(Sf, Pshading, label = 'shading')
plt.plot(Sf, Psheet, label = 'emitter')
plt.legend(loc=0)
plt.show()

df = np.linspace(0.1e-4,2e-4)
Ptotal, Presistivity, Pshading, Psheet  = pv.Pf_total(2.0,0.035, 0.2, 1.6e-6, 50, 100e-4, df, 0.650)
plt.ylim(0,20)
plt.xlabel('finger thickness (µm)')
plt.ylabel('power loss (%)')
plt.plot(df*1e4, Ptotal, label='total')
plt.plot(df*1e4, Presistivity, label='resistive')
#plt.plot(df, Pshading, label = 'shading')
#plt.plot(df, Psheet, label = 'emitter')
plt.legend(loc=0)
plt.show()
