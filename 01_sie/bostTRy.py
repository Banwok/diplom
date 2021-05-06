import Dll1
import matplotlib.pyplot as plt
from math import *
from datetime import datetime
import time
import pickle

#import nolds
#import numpy as np
#from sympy import diff, symbols, cos, sin, exp, simplify

#plt.rcParams["figure.figsize"] = (14,9)
#plt.rcParams.update({'font.size': 40})
#plt.subplots_adjust(left=None, right= None, top=0.98, bottom = 0.15)

a_p = 0.001
a_z = 0.0001
alfa_0 = 0.23
b_p = 0.01
b_z = 0.01
g_p = 0.001
q_p = 6.
q_z = 5.5
Z0 = 0.
Z1 = 1.
P0 = 0.
P1 = 1. 
k_z = 0.15
k_p = 0.05
Q0 = 5. 
Qs = 0.2
t_porog = 100000

t = 200000
t0 = 0
h = 0.1
x0 = 3.7
y0 = 0.2

i_py = 0.0049#
i_end = 0.0077

#print("now ")
#start_time = datetime.now()
#rX = list()
#rY = list()
#rT = list()

#w_q = 0.0063

#Dll1.init(a_p,  a_z,  alfa_0,  b_p,  b_z,  g_p,  q_p,  q_z, Q0,  Qs)
##Dll1.init_rarely(Z0,  Z1,  P0,  P1,  k_z,  k_p, t_porog)
#Dll1.init_runge(x0, y0, t, h, t0)
#Dll1.retuRelease(rX, rY, rT, w_q)

#print("Time C++ : ", datetime.now() - start_time)

#print(type(rX))
#print(len(rX))

#print(type(rY))
#print(len(rY))

#print(type(rT))
#print(len(rT))

##plt.title("Qs = 0.01, θz = 5.5, θp =6, ω = 0.0")
##plt.plot(rT,rX, label = 'График fz(t)')
##plt.plot( rT,rY, label = 'График fp(t)')
#plt.plot( rX,rY)
##plt.xlabel('Ось time ')
##plt.ylabel('Ось Z(t)')
#plt.xlabel('ω$_s$')
#plt.ylabel('Z$_m$$_a$$_x$') 
#plt.legend()

#plt.show()

###################################################################################
#mfu = list()
#tfu = list()
#ly1= list()
#ly2 =list()
#xsh = list()
#with open ('lmax', 'rb') as fp:
#    mfu = pickle.load(fp)
#with open ('retun', 'rb') as fp:
#    tfu = pickle.load(fp)
#with open ('lya1', 'rb') as fp:
#    ly1 = pickle.load(fp)
#with open ('lya2', 'rb') as fp:
#    ly2 = pickle.load(fp)
    
#with open ('t_lya', 'rb') as fp:
#    xsh = pickle.load(fp)

    


#fig = plt.figure(figsize=(10,7))
#ax1 = fig.add_subplot(1,1,1)


## zero line

#ax1.scatter(tfu, mfu, s = 1.5, color= 'darkred')
## ax1.set_xlabel('r')
#ax2 = ax1.twinx()
#ax2.plot(xsh,ly1,label = 'l1')#s = 3, color= 'blue', 
#ax2.plot(xsh,ly2,  label = 'l2')#s = 3, color= 'orange',
##ax1.set_xlim(0.0045, 0.0079)
##ax1.set_ylim(-0.08, 5.5)
#ax1.grid('on')
#plt.show()

########################################
xera = list()
retun = list()

lya1 = list()
lya2 = list()
w_lya = list()

Dll1.init(a_p,  a_z,  alfa_0,  b_p,  b_z,  g_p,  q_p,  q_z, Q0,  Qs)
#Dll1.init_rarely(Z0,  Z1,  P0,  P1,  k_z,  k_p, t_porog)
Dll1.init_runge(x0, y0, t, h, t0)
Dll1.init_wq(i_py, i_end) 

print("Bifurcation begin.. ")
Dll1.bMotherFurkation(xera, retun)
#Dll1.lyapunov_solution(xera, retun, lya1, lya2, w_lya)
#Dll1.seq_Furkation(xera, retun)

print("record to the file")
print("begin")
with open('lmax', 'wb') as fp:
    pickle.dump(xera, fp)

with open('retun', 'wb') as fp:
    pickle.dump(retun, fp)

with open('lya1', 'wb') as fp:
    pickle.dump(lya1, fp)

with open('lya2', 'wb') as fp:
    pickle.dump(lya2, fp)
with open('t_lya', 'wb') as fp:
    pickle.dump(w_lya, fp)

print("end")
sizeL = len(xera)
sizeT = len(retun)

#if(sizeL < sizeT):
#    #retun.pop()
#    print("4ert xera")
#    xera.append(0)
#if(sizeL > sizeT):
#    #xera.pop()
#    print("4ert retun")
#    retun.append(0)

print("LOCAL = ", sizeL)
print("W_Q = ", sizeT)

#print("Qs = 0.01, θz = 5.5, θp =6")

#plt.scatter(retun, xera,  s = 1.5, color= 'darkred')
#plt.scatter(w_lya,lya1,s = 3, color= 'blue', label = 'l1')
#plt.scatter(w_lya,lya2,s = 3, color= 'black', label = 'l2')
##plt.title("Qs = ",Qs, "θz = ", q_z, "θp = ", q_p)
#plt.xlabel('ω$_s$')
#plt.ylabel('Z$_m$$_a$$_x$') 
#plt.show()

fig = plt.figure(figsize=(14,9))
ax1 = fig.add_subplot(1,1,1)
ax1.scatter(retun, xera, s = 1.5, color= 'darkred')
## ax1.set_xlabel('r')
ax2 = ax1.twinx()
ax2.scatter(w_lya,lya1,s = 3, color= 'blue', label = 'l1')#s = 3, color= 'blue', 
ax2.scatter(w_lya,lya2,s = 3, color= 'orange',  label = 'l2')#s = 3, color= 'orange',
##ax1.set_xlim(0.0045, 0.0079)
##ax1.set_ylim(-0.08, 5.5)
ax2.grid('on')
#ax1.xlabel('ω$_s$')
#ax1.ylabel('Z$_m$$_a$$_x$') 
plt.show()

#######################################
#a_p, a_z, alfa_0, b_p, b_z, g_p, q_p, q_z, Z0, Z1, P0, P1, k_z, k_p, Q0, Qs, t, t0, w_q, Z, P = symbols('a_p a_z alfa_0 b_p b_z g_p q_p q_z Z0 Z1 P0 P1 k_z k_p Q0 Qs t t0 w_q Z P')
#fx = (-(a_z + g_p * P) * Z + b_z * (Z0 - (Z0 - Z1) / (1 + exp(-(Q0 + alfa_0 * Z + Qs * sin(w_q * t) - q_z) / k_z))))
#fy = (-a_p * P + b_p * (P0 - (P0 - P1) / (1 + exp(-(Q0 + alfa_0 * Z + Qs * sin(w_q * t) - q_p) / k_p))))

#fx_difZ = diff(fx, Z, 1)
#fx_difP = diff(fx, P, 1)

## simplify(fx_difZ) 
## simplify(fx_difP) 

#print("fx_Z =", fx_difZ)
#print("fx_P =", fx_difP)

#fy_difZ = diff(fy, Z, 1)
#fy_difP = diff(fy, P, 1)

## simplify(fy_difZ) 
## simplify(fy_difP) 

#print("fy_Z = ", fy_difZ)
#print("fy_P = ", fy_difP)
