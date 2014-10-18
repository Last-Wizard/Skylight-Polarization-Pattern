# -*- coding: utf-8 -*-
# 2014.10.16

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from math import pi, sin, cos, acos, atan, sqrt,tan

def rayleigh_2d(Hs=0.0,As=0.0):
	'''
	@Hs: elevation angle of the sun (0,90)
	@As: azimuth angle of the sun (-180,180)
	'''
	# maximum degree of polarization in Rayleigh model
	deltamax = 1
	# normalized radius of celestial sphere
	R = 1
	# transform the azimuth angle from 180 to 360
	# (-pi,-pi/2)->(pi/2,pi),(-pi/2,0)->(pi,pi*3/2)
	# (0,pi/2)->(pi*3/2,2*pi),(pi/2,pi)->(0,pi/2)
	# if As in (pi/2,pi), phis = As-pi/2
	# else phis = As+3/2*pi

	# thetas: zenith angle of the sun (90-Hs)
	thetas = pi/2 - Hs*pi/180.0
	As = As*pi/180
	# phis: azimuth angle of the sun (0,360)
	# relationship between sun and the sphere coordinate
	if pi/2 <= As <= pi:
		phis = As - pi/2
	else:
		phis = As + 1.5*pi

	# dop: degree of polarization
	# aop: angle of polarization	 
	dop = np.zeros((100,400)) # np.zeros(), dtype=float
	aop = np.zeros((100,400))
	x = np.zeros((100,400))
	y = np.zeros((100,400))

	# initial location of the sun
	xs,ys = 0,0
	# for loop	 
	i,j = 0,0

	# Ho: elevation angle of the observer
	# Ao: azimuth angle of the observer
	for Ho in np.linspace(0,pi/2,100,endpoint=True):
		for Ao in np.linspace(0,2*pi,400,endpoint=True):
			thetao = pi/2 - Ho
			phio = Ao # for simplification
			# if pi/2 <= Ao <= pi:
			# 	phio = Ao - pi/2
			# else:
			# 	phio = Ao + 1.5*pi
			# gamma: angular difference between sky point and the sun
			gamma = acos(sin(thetas)*sin(thetao)*cos(phio-phis)+cos(thetas)*cos(thetao))

			# dop: degree of polarization
			dop[i,j] = (deltamax*(sin(gamma))**2)/(1+cos(gamma)**2)

			# aop: angle of polarization
			if sin(phio-phis)*sin(thetas) == 0:
				aop[i,j] = pi/2
			else:
				aop[i,j] = atan((sin(thetao)*cos(thetas)-cos(thetao)*cos(phio-phis)*sin(thetas))/(sin(phio-phis)*sin(thetas)))

			# complex sphere mapping
			x[i,j] = 2*R*tan(thetao/2)*cos(phio)
			y[i,j] = 2*R*tan(thetao/2)*sin(phio)

			j += 1

		i += 1
		j = 0

	# convert angle from radians to degree
	aop = np.rad2deg(aop)

	# location of the sun
	# complex sphere mapping
	xs = 2*R*tan(thetas/2)*cos(phis)
	ys = 2*R*tan(thetas/2)*sin(phis)

	# # # # # # # # # # # # # # # # # # # # #
	# for complex sphere mapping r = 2
	# the equation for this sphere is
	# x^2 + y^2 = 4
	# the equation for (0,0) and (xs,ys) is
	# y = (ys/xs)*xs
	# get the intersection points of the above equations
	meridian_x = np.zeros(2)
	meridian_y = np.zeros(2)

	meridian_x[0] = 2*xs/(sqrt(xs**2 + ys**2))
	meridian_x[1] = - meridian_x[0]

	meridian_y[0] = 2*ys/(sqrt(xs**2 + ys**2))
	meridian_y[1] = -meridian_y[0]


	fig = plt.figure(num='degree of polarizaiton',facecolor='white')
	plt.clf()
	plt.pcolormesh(x,y,dop)
	plt.plot(meridian_x,meridian_y,'w-',xs,ys,'wo',linewidth=2,markersize=8)
	plt.axis('equal')
	plt.axis('off')
	plt.savefig('degree of polarizaiton')
	plt.close()

	fig = plt.figure(num='angle of polarizaiton',facecolor='white')
	plt.clf()
	plt.pcolormesh(x,y,aop)
	plt.plot(meridian_x,meridian_y,'w-',xs,ys,'wo',linewidth=2,markersize=8)
	plt.axis('equal')
	plt.axis('off')
	plt.savefig('angle of polarizaiton')
	plt.close()

# main function
if __name__=='__main__':
	print 'please input the elevation angle of the sun (0,90):',
	Hs = float(raw_input())
	print 'please input the azimuth angle of the sun (-180,180):',
	As = float(raw_input())
	rayleigh_2d(Hs,As)