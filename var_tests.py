#this program calculated the D90, RMS and J STET
#for all stars in a given KELT field and outputs 
#the files, this is for the fields separately
import numpy 
import scipy
from scipy import stats
import math 
from time import time

#import relevant libraries for a list
import glob, os
from os import listdir
from os.path import isfile, join

#which KELT field do you want to reduce?
field = 'field-to-reduce'
#read in the star list
e_name = numpy.loadtxt('path-to-lightcurve-list', dtype = 'string')
w_name = numpy.loadtxt('path-to-lightcurve-list', dtype = 'string')

#get the name holders ready
elc = []
wlc = []

outfile = open('output-file', 'w')
line = 'Name, x, y, RA, Dec, Mean_Mag, RMS, D90, J, K, L, NPTS\n' 
outfile.write(line)
kk = 0

print 'Working on the EAST dataset for '+field+'.'
for ii in range(0, len(e_name)):
	if (kk % 1000 == 0):
		#start the timer
		st = time()

	#get only the variable name
	nm = e_name[ii,2].split('.')
	elc.append(nm[0]+'_mag.data')
	x = e_name[ii,0]
	y = e_name[ii,1]
	ra = e_name[ii,5]
	dec = e_name[ii,6]
	if (os.path.isfile('lightcurve-path') == 1):

		#read in the light curve
		jd, mg, er = numpy.loadtxt('lightcurve-path', unpack = 1)
		
		#do a sigma clip, currently 3 sigma
		clp, lw, upp = scipy.stats.sigmaclip(mg,3, 3)

		#observation number
		nms = len(clp)

		if (nms > 10):	#added in case there are not enough terms
			#get the rms
			rms = numpy.std(clp)

			#get the d90
			cut = long(0.05*nms)
 			clp.sort()	
			d90 = numpy.abs(clp[cut]-clp[-cut])

			#get J stet
			wk = 1.0  #Weighting Factor

			MeanMag = numpy.mean(clp)

			Jt = numpy.arange(nms)*0.0
			Jb = numpy.arange(nms)*0.0
			Kt = numpy.arange(nms)*0.0
			Kb = numpy.arange(nms)*0.0

			for i in range(0, nms-2,2):

				Sigi = (mg[i] - MeanMag)/(er[i])*(numpy.sqrt(nms/(nms-1)))
				Sigj = (mg[i+1] - MeanMag)/(er[i+1])*(numpy.sqrt(nms/(nms-1)))

				Pk = Sigi*Sigj  # pg 853 Stetson 1996 Eq 2 Kinemuchi
				if (Pk > 0.0):
					sgnPk = 1.0
				if (Pk == 0.0): 
					sgnPk = 0.0
				if (Pk < 0.0):
					sgnPk = -1.0

				Jt[i] = wk*sgnPk*(numpy.sqrt(abs(Pk)))	#Kinemuchi eq.1 (Numerator)
				Jb[i] = (wk)				#Kinemuchi eq.1 (Denominator)
				Kt[i] = abs(Sigi)			#Kinemuchi eq.5 (Numerator)
				Kb[i] = abs(Sigi**(2.0))		#Kinemuchi eq.5 (Denominator)


			jstet = sum(Jt)/sum(Jb)        					#Eq 1
			kstet = ((1.0/nms)*sum(Kt))/(numpy.sqrt((1.0/nms)*sum(Kb)))   	#Eq 5
			lstet = jstet*kstet/(0.7908)        				#Eq7

			#print the statistics to a file

			line = str(nm[0]+'_mag.data')+','+str(x)+','+str(y)+','+str(ra)+','+str(dec)+','+str(numpy.mean(clp))+','+str(rms)+','+str(d90)+','+str(jstet)+','+str(kstet)+','+str(lstet)+','+str(nms)+'\n'
			outfile.write(line)
	
			#stop the watch
			kk = kk+1
		
			if (kk % 1000 == 0) and (kk > 0):
				fn = time()
				print 'For field E'+field+', 1000 stars had the stats measured in '+str(fn-st)+'s. '+str(len(e_name)-kk)+' stars remain.'

outfile.close
outfile = open('output-file', 'w')
kk =0
print 'Working on the WEST dataset for '+field+'.'
for ii in range(0, len(w_name)):
	if (kk % 1000 == 0):
		#start the timer
		st = time()

	#get only the variable name
	nm = w_name[ii,2].split('.')
	wlc.append('lightcurve-path')
	x = w_name[ii,0]
	y = w_name[ii,1]
	ra = w_name[ii,5]
	dec = w_name[ii,6]
	if (os.path.isfile('lightcurve-path') == 1):
		
		#read in the light curve
		jd, mg, er = numpy.loadtxt('lightcurve-path', unpack = 1)
		
		#do a sigma clip, currently 3 sigma
		clp, lw, upp = scipy.stats.sigmaclip(mg,3, 3)

		#observation number
		nms = len(clp)

		if (nms > 10): #added in case there are not enough values in the clip
			#get the rms
			rms = numpy.std(clp)

			#get the d90
			cut = long(0.05*nms)
 			clp.sort()	
			d90 = numpy.abs(clp[cut]-clp[-cut])

			#get J stet
			wk = 1.0  #Weighting Factor

			MeanMag = numpy.mean(clp)

			Jt = numpy.arange(nms)*0.0
			Jb = numpy.arange(nms)*0.0
			Kt = numpy.arange(nms)*0.0
			Kb = numpy.arange(nms)*0.0

			for i in range(0, nms-2,2):

				Sigi = (mg[i] - MeanMag)/(er[i])*(numpy.sqrt(nms/(nms-1)))
				Sigj = (mg[i+1] - MeanMag)/(er[i+1])*(numpy.sqrt(nms/(nms-1)))

				Pk = Sigi*Sigj  # pg 853 Stetson 1996 Eq 2 Kinemuchi
				if (Pk > 0.0):
					sgnPk = 1.0
				if (Pk == 0.0): 
					sgnPk = 0.0
				if (Pk < 0.0):
					sgnPk = -1.0

				Jt[i] = wk*sgnPk*(numpy.sqrt(abs(Pk)))	#Kinemuchi eq.1 (Numerator)
				Jb[i] = (wk)				#Kinemuchi eq.1 (Denominator)
				Kt[i] = abs(Sigi)			#Kinemuchi eq.5 (Numerator)
				Kb[i] = abs(Sigi**(2.0))		#Kinemuchi eq.5 (Denominator)


			jstet = sum(Jt)/sum(Jb)        					#Eq 1
			kstet = ((1.0/nms)*sum(Kt))/(numpy.sqrt((1.0/nms)*sum(Kb)))   	#Eq 5
			lstet = jstet*kstet/(0.7908)        				#Eq7

			#print the statistics to a file
			line = str('lightcurve-path')+','+str(x)+','+str(y)+','+str(ra)+','+str(dec)+','+str(numpy.mean(clp))+','+str(rms)+','+str(d90)+','+str(jstet)+','+str(kstet)+','+str(lstet)+','+str(nms)+'\n'
			outfile.write(line)
	
			#stop the watch
			kk = kk+1
		
			if (kk % 1000 == 0) and (kk > 0):
				fn = time()
				print 'For field W'+field+', 1000 stars had the stats measured in '+str(fn-st)+'s. '+str(len(w_name)-kk)+' stars remain.'

outfile.close
