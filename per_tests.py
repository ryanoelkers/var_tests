#this program calculated the D90, RMS and J STET
#for all stars in a given KELT field and outputs 
#the files, this is for the fields separately
import numpy 
import scipy
from scipy import stats
from scipy.signal import spectral
import math 
from time import time
from astropy import stats

#import relevant libraries for a list
import glob, os
from os import listdir
from os.path import isfile, join

#which field do you want to reduce?
field = 'field-to-reduce'

#read in the light curve names for the EAST
evars = numpy.loadtxt('path-to-lightcurve-list', skiprows =1 ,dtype = 'string')

#find the light curve with the most points
for ii in range(0, len(evars)):
	lst = evars[ii].split(',')
	npts = float(lst[11])
	if (ii == 0):
		mxpts = npts
		idx = 0
	else:
		if (ii > 0):
			if (npts > mxpts):
				mxpts = npts
				idx = ii
#now read in the longest light curve and get the period
nm = evars[idx].split(',')
nm1 = nm[0].split('_')
jd, mag, er = numpy.loadtxt('path-to-lightcurve', unpack = 1)

#min and max peaks
mxp = long((jd[-1]-jd[0]))
mnp = 0.1

#read in the star list
e_name = numpy.loadtxt('path-to-lightcurve-list', dtype = 'string')
w_name = numpy.loadtxt('path-to-lightcurve-list', dtype = 'string')

#get the name holders ready
elc = []
wlc = []

#open the output and write the file header
outfile = open('out-put-file', 'w')
line = 'Name, x, y, RA, Dec, Mean_Mag, Std_Mag, P1, S1, P2, S2, P3, S3, P4, S4,P5, S5,P6, S6,P7, S7,P8, S8, P9, S9, P10, S10\n' 
outfile.write(line)
kk = 0

print 'Working on the EAST dataset for '+field+'.'
for jj in range(0, len(e_name)):
	#start or restart the timer based on your current iteration
	if (kk % 1000 == 0):
		#start the timer
		st = time()

	#get only the variable name
	nm = e_name[jj,2].split('.')
	nm1 = nm[0].split('_')
	elc.append('light-curve-name')
	x = e_name[jj,0]
	y = e_name[jj,1]
	ra = e_name[jj,5]
	dec = e_name[jj,6]
	#print elc[-1]
	#check to make sure the file actually exists
	if (os.path.isfile('light-curve-path') == 1):
		
		#read in the light curve
		jd, mg, er = numpy.loadtxt('light-curve-path', unpack = 1)
		#ers = 
		#do Lomb-Scargle
		freq, pwr = stats.LombScargle(jd, mg).autopower(minimum_frequency=1./mxp, maximum_frequency=1./mnp)

		pks = []
		pkpw = []
		#get the 10 largest peaks
		for ii in range(1, len(pwr)-1):
			if (pwr[ii] > pwr[ii-1]) and (pwr[ii] > pwr[ii+1]):
				pks.append(freq[ii])
				pkpw.append(pwr[ii])
		#sort the peaks based on power
		tt = zip(*sorted(zip(pkpw, pks)))[0]
		vv = zip(*sorted(zip(pkpw,pks)))[1]

		#write the line with the information
		line = str('light-curve-path')+','+str(x)+','+str(y)+','+str(ra)+','+str(dec)+','+str(numpy.mean(mg))+','+str(numpy.std(mg))+','+str(1.0/vv[-1])+','+str(tt[-1])+','+str(1.0/vv[-2])+','+str(tt[-2])+','+str(1.0/vv[-3])+','+str(tt[-3])+','+str(1.0/vv[-4])+','+str(tt[-4])+','+str(1.0/vv[-5])+','+str(tt[-5])+','+str(1.0/vv[-6])+','+str(tt[-6])+','+str(1.0/vv[-7])+','+str(tt[-7])+','+str(1.0/vv[-8])+','+str(tt[-8])+','+str(1.0/vv[-9])+','+str(tt[-9])+','+str(1.0/vv[-10])+','+str(tt[-10])+'\n'
		outfile.write(line)

	kk = kk+1
	if (kk % 1000 == 0):
		ed = time()
		print 'Done with 1000 in '+str(ed-st)+'. Working on the next 1000. '+str(len(e_name)-kk)+' stars remain in E'+field+'.'
#close the file
outfile.close

#read in the light curve names for the WEST
wvars = numpy.loadtxt('path-to-lightcurve-list', skiprows=1,dtype = 'string')

#find the light curve with the most points
for ii in range(0, len(wvars)):
	lst = wvars[ii].split(',')
	npts = float(lst[11])
	if (ii == 0):
		mxpts = npts
		idx = 0
	else:
		if (ii > 0):
			if (npts > mxpts):
				mxpts = npts
				idx = ii
#read in the longest light curve
nm = wvars[idx].split(',')
nm1 =nm[0].split('_')
jd, mag, er = numpy.loadtxt('path-to-lightcurve-list', unpack = 1)
#min and max peaks
mxp = long((jd[-1]-jd[0]))
mnp = 0.1

#open the output and write the file header
outfile = open('path-to-lightcurve-list', 'w')
line = 'Name, x, y, RA, Dec, Mean_Mag, Std_Mag, P1, S1, P2, S2, P3, S3, P4, S4,P5, S5,P6, S6,P7, S7,P8, S8, P9, S9, P10, S10\n' 
outfile.write(line)
kk = 0

print 'Working on the WEST dataset for '+field+'.'
for jj in range(0, len(w_name)):
	#start or restart the timer based on your current iteration
	if (kk % 1000 == 0):
		#start the timer
		st = time()

	#get only the variable name
	nm = w_name[jj,2].split('.')
	nm1 = nm[0].split('_')
	wlc.append('light-curve-path')
	x = w_name[jj,0]
	y = w_name[jj,1]
	ra = w_name[jj,5]
	dec = w_name[jj,6]
		
	#check to make sure the file actually exists
	if (os.path.isfile('light-curve-path') == 1):
		
		#read in the light curve
		jd, mg, er = numpy.loadtxt('light-curve-path', unpack = 1)
		
		#do Lomb-Scargle
		freq, pwr = stats.LombScargle(jd, mg).autopower(minimum_frequency=1./mxp, maximum_frequency=1./mnp)

		pks = []
		pkpw = []
		#get the 10 largest peaks
		for ii in range(1, len(pwr)-1):
			if (pwr[ii] > pwr[ii-1]) and (pwr[ii] > pwr[ii+1]):
				pks.append(freq[ii])
				pkpw.append(pwr[ii])
		#sort the peaks based on power
		tt = zip(*sorted(zip(pkpw, pks)))[0]
		vv = zip(*sorted(zip(pkpw,pks)))[1]

		#write the line with the information
		line = str('light-curve-path')+','+str(x)+','+str(y)+','+str(ra)+','+str(dec)+','+str(numpy.mean(mg))+','+str(numpy.std(mg))+','+str(1.0/vv[-1])+','+str(tt[-1])+','+str(1.0/vv[-2])+','+str(tt[-2])+','+str(1.0/vv[-3])+','+str(tt[-3])+','+str(1.0/vv[-4])+','+str(tt[-4])+','+str(1.0/vv[-5])+','+str(tt[-5])+','+str(1.0/vv[-6])+','+str(tt[-6])+','+str(1.0/vv[-7])+','+str(tt[-7])+','+str(1.0/vv[-8])+','+str(tt[-8])+','+str(1.0/vv[-9])+','+str(tt[-9])+','+str(1.0/vv[-10])+','+str(tt[-10])+'\n'
		outfile.write(line)

	kk = kk+1
	if (kk % 1000 == 0):
		ed = time()
		print 'Done with 1000 in '+str(ed-st)+'. Working on the next 1000. '+str(len(w_name)-kk)+' stars remain in W'+field+'.'
#close the file
outfile.close
