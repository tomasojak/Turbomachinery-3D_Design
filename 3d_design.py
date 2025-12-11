import numpy as np 
import matplotlib.pyplot as plt

#given values 
mass_flow  =  35
rpm  =  16000
U_tip = 430
r_h_r_t = 0.375
rho = 1.204
kappa = 1.4
R_gas = 287
T = 273+3
R_midspan = 0.3 #stand in value for now
Phi_midspan = 0.5


omega = 2*np.pi/60*rpm
r_tip = U_tip/omega
r_hub = r_tip*r_h_r_t
U_hub  =  omega*r_hub
print(f'r hub {r_hub} m')
print(f'r tip {r_tip} m')

r_midspan = (r_tip+r_hub) / 2

c_ax = mass_flow / np.pi / (np.power( r_tip,2 ) - np.power( r_hub,2 )) / rho
print(c_ax)

#relative Mach at inlet - need to change c_ax to w 
#this is absolute Mach at inlet - relative is calculated lower 
mach = c_ax/np.sqrt(kappa*T*R_gas)

print(f'Mach number at inlet {mach}')

#calculating 
alfa_1 = 0
# beta_1_tip = np.arctan(U_tip/c_ax)
# beta_1_hub = np.arctan(U_hub/c_ax)

#A2 = r_midspan*((0.5-R_midspan)*2*omega*r_midspan+omega*r_midspan)
A2  =  r_midspan * ((-2*R_midspan + 1) * omega * r_midspan)


def three_d_design(A2, r, i):
    U = r*omega
    c_theta_2 = A2/r
    beta_1 = np.arctan(U/c_ax)
    
    #beta_2 = np.arctan((0.5-R)*2*r*omega/c_ax)
    beta_2 = np.arctan((-c_theta_2 + U) / c_ax)
    #alfa_2 = np.arctan((np.tan(beta_2)*c_ax+r*omega)/c_ax)
    #alfa_2 = np.arctan(omega*r/c_ax-np.tan(beta_2))
    alfa_2  =  np.arctan((U - c_ax*np.tan(beta_2)) / c_ax)
    deg = 180/np.pi
    beta_1_array[i] = beta_1*deg
    beta_2_array[i] = beta_2*deg
    alfa_2_array[i] = alfa_2*deg
    radius[i] = r

delta = 0.01
nr_iter = int(np.round((r_tip-r_hub)/delta))
print(nr_iter)

beta_1_array = np.zeros(nr_iter)
beta_2_array = np.zeros(nr_iter)
alfa_1_array = np.zeros(nr_iter)
alfa_2_array = np.zeros(nr_iter)
radius = np.zeros(nr_iter)
Mach_array = np.zeros(nr_iter)

for i in range(nr_iter):
    three_d_design(A2,r_hub+i*delta, i)
    alfa_1_array[i] = alfa_1
    w = np.sqrt(np.power(r_hub+i*delta*omega,2)+np.power(c_ax,2))
    Mach_array[i] = w/np.sqrt(kappa*T*R_gas)

print(Mach_array)
##### PLOTTING #####

fontSize  =  10

plt.rcParams.update({
"font.family": "sans-serif",
"font.size": fontSize,
"axes.labelsize": fontSize,
"axes.titlesize": fontSize,
"legend.fontsize": fontSize,
"xtick.labelsize": fontSize,
"ytick.labelsize": fontSize,
})


# plot
fig, ax  =  plt.subplots()

ax.plot(radius, alfa_1_array, label = 'alfa 1')
ax.plot(radius, alfa_2_array, label = 'alfa 2')
ax.plot(radius, beta_1_array, label = 'beta 1')
ax.plot(radius, beta_2_array, label = 'beta 2')

# ax.set(xlim = (0, 8), xticks = np.arange(1, 8),
#        ylim = (0, 8), yticks = np.arange(1, 8))
plt.legend()
plt.show()