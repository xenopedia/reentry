
import math
import matplotlib.pyplot as plt
import numpy as np

import usatmos
import rk

GS = 0.295912208285591095E-03
GE = 3.98600433E+14
ERAD = 6378136.6

###Apollo###
Cd = 1.35 #aerodynamic coeff
A = 12 #cross-sectional area
mc = 5300.0 #mass
LoverD = 0.268
###########


###Hayabusa###
#Cd = 1.3 #aerodynamic coeff
#A = 0.128 #cross-sectional area
#mc = 18.0 #mass
#LoverD = 0.0
#############


def fall_force_model (t,y,eq):

    derivs = [None]*4
    gd = math.sqrt ( y[0]*y[0] + y[2]*y[2]  );
    gv = math.sqrt ( y[1]*y[1] + y[3]*y[3]  );
    gd3 = gd*gd*gd

    derivs[0] = y[1];
    derivs[2] = y[3];

    density, pressure, temperature = usatmos.Atmosphere((gd-ERAD)/1000.0)

    derivs[1] = -GE/gd3*y[0] - 0.5*density*gv*y[1]*Cd*A/mc
    derivs[3] = -GE/gd3*y[2] - 0.5*density*gv*y[3]*Cd*A/mc + 0.5*density*gv*y[1]*LoverD*Cd*A/mc

    return derivs[eq]

speed = 11040.0
angle = 6.52
pv = [0, speed*math.cos(angle*math.pi/180.), 121920+ERAD, -speed*math.sin(angle*math.pi/180.)] #apollo ENTRY interface
time = 0.0
step = 0.1

T = []
X = []
Y = []
VX = []
VY = []
G = []
H = []
LD = []

#earth surface
XE = []
YE = []
DOWNRANGE = []

#control variables
iv = [None]*4 #intermediary run
integral = 0 
pe = 0 #prev error
nrseq = 0

while math.sqrt(pv[0]*pv[0]+pv[2]*pv[2]) > ERAD:
    newstep = rk.rk45(fall_force_model, time, pv, step, 1E-7)

    time += step

    T.append(time)
    X.append(pv[0]/1000.)
    VX.append(pv[1])
    Y.append(pv[2]/1000.)
    VY.append(pv[3])
    H.append((math.sqrt(pv[0]*pv[0]+pv[2]*pv[2])-ERAD)/1000.0)
    LD.append(LoverD)

    if time > step*3:
        G.append( (math.sqrt(VX[-1]*VX[-1] + VY[-1]*VY[-1]) - math.sqrt(VX[-2]*VX[-2] + VY[-2]*VY[-2]) )/step/9.81  )
    else:
        G.append (0)
    
    step = newstep

    angle = math.pi/2 - math.atan2(pv[2], pv[0])
    if angle<0:
       angle = math.pi/2 - angle
    downrange_dist = ERAD * angle / 1000.

    XE.append(ERAD*math.cos(math.atan2(pv[2], pv[0]))/1000.)
    YE.append(ERAD*math.sin(math.atan2(pv[2], pv[0]))/1000.)
    DOWNRANGE.append(downrange_dist)


# H vs. T
#    plt.plot(T, H, 'b', lw=3, alpha=0.8)
#    plt.grid(True, which='major')
#    plt.xlim(T[-1]-10,T[-1]+10)
#    plt.ylim(max(H[-1]-10,0), H[-1]+10)
#    plt.plot(T[-1], H[-1], 'ro', lw=2, alpha=1.)
#    name = "img" + str(10000+nrseq) + ".png"
#    plt.savefig(name)
#    plt.clf()


# Height vs. Downrange
#    plt.plot(DOWNRANGE, H, 'b', lw=4, alpha=0.6)
#    plt.grid(True, which='major')
#    plt.xlim(DOWNRANGE[-1]-30,DOWNRANGE[-1]+30)
#    plt.ylim(max(H[-1]-10,0), H[-1]+10)
#    plt.plot(XE, YE, 'black', lw=3, alpha=0.8)
#    plt.plot(DOWNRANGE[-1], H[-1], 'ro', lw=2, alpha=1., ms=10)
#    plt.xlabel("Downrange (km)")
#    plt.ylabel("Height (km)")
#    plt.title("Time: "+str(int(time))+" (s)")
#    name = "img" + str(100000+nrseq) + ".png"
#    plt.savefig(name)
#    plt.clf()

final_angle = math.pi/2 - math.atan2(pv[2], pv[0])
if final_angle<0:
    final_angle = math.pi/2 - final_angle
final_dist = ERAD * final_angle

print ("final distance ", final_dist/1000.0)

plt.figure(figsize=(8,6))
plt.subplot(2,2,1)
plt.plot(T, H, 'b', lw=3, alpha=0.8)
plt.grid(True, which='major')
plt.xlabel('Time (s)')
plt.ylabel('Height (km)')
plt.title(' ')

plt.subplot(2,2,2)
plt.plot(T, G, 'r', lw=3, alpha=0.8)
#plt.axis([0,200,0,200])
plt.grid(True, which='major')
plt.xlabel('Time (s)')
plt.ylabel('g')
plt.title(' ')

plt.subplot(2,2,3)
plt.plot(DOWNRANGE, H, 'g', lw=3, alpha=0.8)
#plt.plot(XE, YE, 'black', lw=3, alpha=0.8)
#plt.axis([0,200,0,200])
plt.grid(True, which='major')
plt.xlabel('Range (km)')
plt.ylabel('Height (km)')
plt.title(' ')

plt.subplot(2,2,4)
plt.plot(T, LD, 'black', lw=3, alpha=0.8)
#plt.axis([0,200,0,200])
plt.grid(True, which='major')
plt.xlabel('Time (s)')
plt.ylabel('L/D')
plt.title(' ')

##plt.xticks(np.arange(0,200,10))
##plt.yticks(np.arange(0,200,10))
##plt.minorticks_on()
plt.tight_layout()
plt.savefig("A8.svg")
#plt.show()
