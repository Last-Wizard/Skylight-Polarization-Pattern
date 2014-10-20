# -*- coding: utf-8 -*-
# wang bo (wangbo.hfut@outlook.com)
# hfut
# 2014.03.27
# reference: 
# http://www.oceanopticsbook.info/view/radiative_transfer_theory/level_2/the_vector_radiative_transfer_equation

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from math import *
import numpy as np

# maximum degree of polarization
deltamax = 1

# normalized radius of celestial sphere
R = 1

# elevation angle of the sun
Hs = pi/6

# azimuth angle of the sun
As = 3*pi/5
# (-pi,-pi/2)->(pi/2,pi),(-pi/2,0)->(pi,pi*3/2)
# (0,pi/2)->(pi*3/2,2*pi),(pi/2,pi)->(0,pi/2)
# if x in (pi/2,pi), y = x-pi/2
# else y = x+3/2*pi

thetas = pi/2 - Hs

# relationship between sun and the sphere coordinate
if pi/2 <= As <= pi:
	phis = As - pi/2
else:
	phis = As + 1.5*pi

# coordinate of observer
# Ho: elevation angle of observer
# thetao = pi/2 - Ho
# Ao: azimuth angle of observer

# degree of polarization
dop = np.zeros((100,400))
# angle of polarization
aop = np.zeros((100,400))

x = np.zeros((100,400))
y = np.zeros((100,400))
z = np.zeros((100,400))

# Stokes vector of natural light
S_in = np.array([1,0,0,0],dtype=float)

# 4 different Stokes vector of scattering light
# S_out = [I,Q,U,V]^T
I = np.zeros((100,400))
Q = np.zeros((100,400))
U = np.zeros((100,400))
V = np.zeros((100,400))

# rotation matrix
R_in = np.zeros((4,4))
R_out = np.zeros((4,4))

# reduced muller matrix for Rayleigh scattering
M = np.zeros((4,4))

# phase matrix
# P = R_out * M * R_in
P = np.zeros((4,4))

i = 0
j = 0

for Ho in np.linspace(0,pi/2,100,endpoint=False):
	for Ao in np.linspace(0,2*pi,400,endpoint=True):
		thetao = pi/2 - Ho
		phio = Ao

#		if pi/2 <= Ao <= pi:
#			phio = Ao - pi/2
#		else:
#			phio = Ao + 1.5*pi
		
		# gamma: angular difference between sky element and sun
		gamma = acos(sin(thetas)*sin(thetao)*cos(phio-phis)+cos(thetas)*cos(thetao))

		if 0 <= phis <= pi:
			if 0 <= (Ao-phis) <= pi:
				alpha_in = acos((-cos(thetao)+cos(thetas)*cos(gamma))/(sin(gamma)*sin(thetas)))
				alpha_out = acos((-cos(thetas)+cos(thetao)*cos(gamma))/(sin(gamma)*sin(thetao)))
			else:
				alpha_in = acos((-cos(thetao)+cos(thetas)*cos(gamma))/(-sin(gamma)*sin(thetas)))
				alpha_out = acos((-cos(thetas)+cos(thetao)*cos(gamma))/(-sin(gamma)*sin(thetao)))
		else:
			if (0 <= (Ao-phis) <= pi) or ( 0 <= Ao <= (phis+pi)%(2*pi)):
				alpha_in = acos((-cos(thetao)+cos(thetas)*cos(gamma))/(sin(gamma)*sin(thetas)))
				alpha_out = acos((-cos(thetas)+cos(thetao)*cos(gamma))/(sin(gamma)*sin(thetao)))
			else:
				alpha_in = acos((-cos(thetao)+cos(thetas)*cos(gamma))/(-sin(gamma)*sin(thetas)))
				alpha_out = acos((-cos(thetas)+cos(thetao)*cos(gamma))/(-sin(gamma)*sin(thetao)))

		if gamma < 10**-4:
			alpha_out = 0
		else:
			pass

		if thetao < 10**-4:
			if 0 <= phis <= pi:
				if 0 <= (Ao-phis) <= pi:
					alpha_out = Ao-phis
				else:
					alpha_out = phis-Ao
			else:
				if (0 <= (Ao-As) <= pi) or ( 0 <= Ao <= (phis+pi)%(2*pi)):
					alpha_out = Ao-phis
				else:
					alpha_out = phis-Ao
		else:
			pass

		# alpha_in: angle between (h,v,r) and (s,p,r)
		# alpha_in = asin(sin(thetao)*sin(phis-phio)/sin(gamma))
		# alpha_out: angle between (h',v',r') and (s',p',r')
		# alpha_out = asin(sin(thetas)*sin(phis-phio)/sin(gamma))

		# rotation matrix
		R_out[0,:] = [1,0,0,0]
		R_out[1,:] = [0,cos(2*alpha_out),-sin(2*alpha_out),0]
		R_out[2,:] = [0,sin(2*alpha_out),cos(2*alpha_out),0]
		R_out[3,:] = [0,0,0,1]

		R_in[0,:] = [1,0,0,0]
		R_in[1,:] = [0,cos(2*alpha_in),-sin(2*alpha_in),0]
		R_in[2,:] = [0,sin(2*alpha_in),cos(2*alpha_in),0]
		R_in[3,:] = [0,0,0,1]

		# reduced muller matrix for Rayleigh scattering
		M[0,:] = [1,-sin(gamma)**2/(1+cos(gamma)**2),0,0]
		M[1,:] = [-sin(gamma)**2/(1+cos(gamma)**2),1,0,0]
		M[2,:] = [0,0,2*cos(gamma)/(1+cos(gamma)**2),0]
		M[3,:] = [0,0,0,2*cos(gamma)/(1+cos(gamma)**2)]

		# phase matrix
		# P = R_out * M * R_in
		R_out = np.dot(R_out,M)
		P = np.dot(R_out,R_in)

		# 4 different Stokes vector of scattering light
		# S_out = [I,Q,U,V]^T
		I[i,j] = np.dot(P[0,:],S_in)
		Q[i,j] = np.dot(P[1,:],S_in)
		U[i,j] = np.dot(P[2,:],S_in)
		V[i,j] = np.dot(P[3,:],S_in)

		aop[i,j] = 0.5*atan(U[i,j]/Q[i,j])

		if Q[i,j] < 0:
			if U[i,j] > 0:
				aop[i,j] += pi/2
			else:
				aop[i,j] -= pi/2
		else:
			pass

		dop[i,j] = sqrt(Q[i,j]**2+U[i,j]**2+V[i,j]**2)/I[i,j]

		x[i,j] = R*sin(thetao)*cos(phio)
		y[i,j] = R*sin(thetao)*sin(phio)

		j += 1

	i += 1
	j = 0

# location of the sun
xs=R*sin(thetas)*cos(phis)
ys=R*sin(thetas)*sin(phis)

# # # # # # # # # # # # # # # # # # # # #
# for directly sphere mapping r = 1
# the equation for this sphere is
# x^2 + y^2 = 1
# the equation for (0,0) and (xs,ys) is
# y = (ys/xs)*xs

# get the intersection points of the above equations
meridian_x = np.zeros(2)
meridian_y = np.zeros(2)

meridian_x[0] = xs/(sqrt(xs**2 + ys**2))
meridian_x[1] = - meridian_x[0]

meridian_y[0] = ys/(sqrt(xs**2 + ys**2))
meridian_y[1] = -meridian_y[0]


font = {'family' : 'Times New Roman',
		'size'   : '14'}

matplotlib.rc('font',**font)

# I = 1
# fig = plt.figure(num='distribution of I',facecolor='white')
# plt.clf()
# plt.pcolormesh(x,y,I)
# plt.plot(meridian_x,meridian_y,'w-',xs,ys,'wo',linewidth=2,markersize=8)
# plt.axis('equal')
# plt.colorbar()

# Q
fig = plt.figure(num='distribution of Q',facecolor='white')
plt.clf()
plt.pcolormesh(x,y,Q,cmap=cm.rainbow)
plt.plot(meridian_x,meridian_y,'w-',xs,ys,'wo',linewidth=2,markersize=8)
plt.axis('equal')
plt.colorbar()

# U
fig = plt.figure(num='distribution of U',facecolor='white')
plt.clf()
plt.pcolormesh(x,y,U,cmap=cm.rainbow)
plt.plot(meridian_x,meridian_y,'w-',xs,ys,'wo',linewidth=2,markersize=8)
plt.axis('equal')
plt.colorbar()

# V = 0
# fig = plt.figure(num='distribution of V',facecolor='white')
# plt.clf()
# plt.pcolormesh(x,y,V)
# plt.plot(meridian_x,meridian_y,'w-',xs,ys,'wo',linewidth=2,markersize=8)
# plt.axis('equal')
# plt.colorbar()

# degree of polarization
fig = plt.figure(num='degree of polarization',facecolor='white')
plt.clf()
plt.pcolormesh(x,y,dop,cmap=cm.jet)
plt.plot(meridian_x,meridian_y,'w-',xs,ys,'wo',linewidth=2,markersize=8)
plt.axis('equal')
plt.colorbar()

# convert angle from radians to degree
aop = np.rad2deg(aop)
# angle of polarizaiton
fig = plt.figure(num='angle of polarization',facecolor='white')
plt.clf()
plt.pcolormesh(x,y,aop,cmap=cm.jet)
plt.plot(meridian_x,meridian_y,'w-',xs,ys,'wo',linewidth=2,markersize=8)
plt.axis('equal')
plt.colorbar()

plt.show()
