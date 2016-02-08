#!/usr/bin/env python
#Program to obtain snrs in frequency regime and of cross correlations using noise windows
#Version: 28/01/2016
#Importing the modules used within this program
from obspy.core import read
import os
import csv
from numpy import correlate
from numpy import argmax
from numpy import corrcoef
import numpy as np
import math
from collections import deque
import matplotlib.pyplot as plt
cdir = os.getcwd()

#prepwork appends zeroes, rotates, takes the fourier transform and extracts the power spectrum
#as well as the conjugate	
def prepwork(corrvals,tw):
	corrvals.extend([0]*(10000-len(corrvals)))
	corrvalsmov = deque(corrvals)
	corrvalsmov.rotate(-tw)
	corrvalsfft = np.fft.fft(corrvalsmov)
	corrvalspow = np.abs(corrvalsfft)
	conjugate = np.conj(corrvalsfft)
	return corrvalsfft, conjugate, corrvalspow

#This function trims a trace (in this case a copy of the original trorig), tapers it, and outputs the trace as a numpy array
#As well as outputting the array as a list (slepwin) for use in prepwork
def sleptap(tror, lowtime, hitime, slepper, slepwid):
	start = tror[0].stats.starttime
	sleptr = tror.copy().trim(start+(lowtime*0.01), start+(hitime*0.01))
	sleptr.taper(max_percentage=slepper,type='slepian',max_length=None,side='both',width=slepwid)
	sleparr = sleptr[0]
	slepwin = sleparr[:].tolist()
	return sleparr, slepwin
	
#Setting the maximum frequency to use in the filter on normalised CC stage
mfreq = 10
maxfreq = str(mfreq)

#Setting empty lists that will be used to store information
lilist = []
biglist = []
stats = []
wholist = []
listoff = []
nmaxlist = []
smax = []
#Read in station list from file set as below
with open('stlist.csv', 'rb') as f:
	reader = csv.reader(f)
	for row in reader:
		stats.append(row)
		
#With station list you can cycle through and set up lists of files which can be 
#cross correlated
#Getting current directory
cdir2 = os.getcwd()
#Looping through and creating two separate lists, one of smaller quake traces, and one with bigger
for i in range(0, len(stats)):
	temp = str(stats[i][0])
	sacfiles = [ f for f in os.listdir(cdir2) if temp in f if f.endswith('.filt')]
	if len(sacfiles) > 1:
		biglist.append(sacfiles[0])
		lilist.append(sacfiles[1])

#More empty lists that will be used
sigunnorm = []
nois1unnorm = []
nois2unnorm = []
nois3unnorm = []
allsnrslil = []
allsnrsbig = []
sumsigcoh = []
sumnois1coh = []
sumnois2coh = []
sumnois3coh = []
s1n1 = []
s2n2 = []
s1s2n1n2 = []

#Slep width set here
slepwid = 0.025
slepper = 0.05
#As we will be making multiple figures in a loop, we need a counter
fignum = 1
#Reading in traces and setting phase lags and lowtimes/hightimes, also performing unnormalised cross correlations
#And pulling out phase coherence values
for j in range(0, len(biglist)):
	stime = -100
	#Copy traces so you don't overwrite them
	trorig = read(lilist[j])
	tr2orig = read(biglist[j])
	listoff.append((lilist[j], biglist[j]))
	tr = trorig.copy()
	tr2 = tr2orig.copy()

	#Applying processing to trace no. 2, this consists of a simple taper, filter and cut
	arr1 = tr[0]
	arr2 = tr2[0]
	t71 = arr1.stats.sac.t7 
	t72 = arr2.stats.sac.t7
	low2 = t72
	temp = tr2[0].stats.starttime
	hi2 = t72 +200
	tr2cop = tr2.copy()
	tr2cop.trim(temp+(low2*0.01),temp+(hi2*0.01))
	tr2cop.filter("bandpass", freqmin=1, freqmax=mfreq, zerophase=True)
	tr2cop.taper(0.05,'cosine',None,'right')
	ptarr2 = tr2cop[0]
	corrlist = []
	wlist = []
	stlist = []
	#Perform normalised cross correlation, I will use the maximum of this to align my unnormalised cross correlations
	for k in range(0,201):
		#Applying processing to trace no. 1
		trloop = tr.copy()
		low1 = t71 + stime
		hi1 = t71 + 200 + stime
		temp1 = tr[0].stats.starttime
		#Processing trace
		trloop.trim(temp1+(low1*0.01),temp1+(hi1*0.01))
		trloop.filter("bandpass", freqmin=1, freqmax=mfreq, zerophase=True)
		trloop.taper(0.05,'cosine',None,'right')
		ptarr1 = trloop[0]
		
		#The normalised cross correlation itself, pulling out coefficient values as a function of time shift (stime)
		fullcorr = corrcoef(ptarr1, ptarr2, rowvar=1)
		corrlist.append(float(fullcorr[1][0]))
		stlist.append(stime)
		stime += 1
		trloop = []
	maxnorm = argmax(corrlist) # This is the maximum of the normalised cross correlation, location in the array
	nmax = max(corrlist)
	nmaxlist.append(nmax)
	#Thus giving us the best lag time
	maxstime = stlist[maxnorm]
	smax.append(maxstime)
	#SO now we run the unnormalised cross correlation, which requires mode to be set to full
	#First changing our time windows
	slow1 = t71 + maxstime
	shi1 = t71 + maxstime+800
	low2 = t72
	hi2 = t72+800
	
	#Need to trim the traces and taper them before correlating, using a DPSS this time
	taparr1, sigwindow = sleptap(trorig, slow1, shi1, slepper, slepwid)
	taparr2, sigwindow2 = sleptap(tr2orig, low2, hi2, slepper, slepwid)

	sigunnorm.append(correlate(taparr1, taparr2, mode='full')) 
	#And we also want to run unnormalised cross correlations on noise windows, placing each one in a list of lists
	nois1low1 = t71 + maxstime - 2700
	nois1hi1 = t71 + maxstime - 1900
	nois1low2 = t72 - 2700
	nois1hi2 = t72 - 1900
	nois2low1 = t71 + maxstime - 1800
	nois2hi1 = t71 + maxstime - 1000
	nois2low2 = t72 - 1800
	nois2hi2 = t72 - 1000
	nois3low1 = t71 + maxstime - 900
	nois3hi1 = t71 + maxstime - 100
	nois3low2 = t72 - 900
	nois3hi2 = t72 - 100
	
	#Performing trimming and tapering on noise windows
	nois1arr1, nois1win = sleptap(trorig, nois1low1, nois1hi1, slepper, slepwid)
	nois1arr2, nois1win2 = sleptap(tr2orig, nois1low2, nois1hi2, slepper, slepwid)
	nois2arr1, nois2win = sleptap(trorig, nois2low1, nois2hi1, slepper, slepwid)
	nois2arr2, nois2win2 = sleptap(tr2orig, nois2low2, nois2hi2, slepper, slepwid)
	nois3arr1, nois3win = sleptap(trorig, nois3low1, nois3hi1, slepper, slepwid)
	nois3arr2, nois3win2 = sleptap(tr2orig, nois3low2, nois3hi2, slepper, slepwid)

	#Performing correlations on noise windows
	nois1unnorm.append(correlate(nois1arr1, nois1arr2, mode = 'full'))
	nois2unnorm.append(correlate(nois2arr1, nois2arr2, mode = 'full'))
	nois3unnorm.append(correlate(nois3arr1, nois3arr2, mode = 'full'))

	
	#One final step we want to perform at this stage is to find the signal to noise ratios in the frequency domain for each trace	
	#Taking the spectra of these windows
	sigwinfft, sigwinconj, sigwinpow = prepwork(sigwindow,800)
	nois1winfft, nois1winconj, nois1winpow = prepwork(nois1win,800)
	nois2winfft, nois2winconj, nois2winpow = prepwork(nois2win,800)
	nois3winfft, nois3winconj, nois3winpow = prepwork(nois3win,800)
	sigwinfft2, sigwinconj2, sigwinpow2 = prepwork(sigwindow2,800)
	nois1winfft2, nois1winconj2, nois1winpow2 = prepwork(nois1win2,800)
	nois2winfft2, nois2winconj2, nois2winpow2 = prepwork(nois2win2,800)
	nois3winfft2, nois3winconj2, nois3winpow2 = prepwork(nois3win2,800)
	
	#Calculating the phase coherence between each pair of stations for the signal and noise windows
	sigcoh = []
	nois1coh = []
	nois2coh = []
	nois3coh = []
	for cn in range(0,len(sigwinfft)):
		sigcoh.append(np.real((sigwinfft[cn]*sigwinconj2[cn])/abs(sigwinfft[cn]*sigwinconj2[cn])))
		nois1coh.append(np.real((nois1winfft[cn]*nois1winconj2[cn])/abs(nois1winfft[cn]*nois1winconj2[cn])))
		nois2coh.append(np.real((nois2winfft[cn]*nois2winconj2[cn])/abs(nois2winfft[cn]*nois2winconj2[cn])))
		nois3coh.append(np.real((nois3winfft[cn]*nois3winconj2[cn])/abs(nois3winfft[cn]*nois3winconj2[cn])))
	
	#Averaging noise window values for use in creating the SNRs (or in this case, NSRs)
	aveofnoiswin = [((a+b+c)/3) for a,b,c in zip(nois1winpow, nois2winpow, nois3winpow)]
	aveofnoiswin2 = [((a+b+c)/3) for a,b,c in zip(nois1winpow2, nois2winpow2, nois3winpow2)]
	
	#Generating frequencies and sampling information every 0.125Hz, also finding nsrs
	fc = 0
	num = 10000
	timestep = 0.01
	freqs = np.fft.fftfreq(num, d = timestep) #This generates the frequency range 
	frang = []
	snrtracelist = []
	snrtracelist2 = []
	multsnrs = []
	sigcohrang = []
	nois1cohrang = []
	nois2cohrang = []
	nois3cohrang = []
	while fc in range(0,3000):
		frang.append(freqs[fc])
		snrtracelist.append(aveofnoiswin[fc]/sigwinpow[fc])
		snrtracelist2.append(aveofnoiswin2[fc]/sigwinpow2[fc])
		multsnrs.append((aveofnoiswin[fc]/sigwinpow[fc])*(aveofnoiswin2[fc]/sigwinpow2[fc]))
		sigcohrang.append(sigcoh[fc])
		nois1cohrang.append(nois1coh[fc])
		nois2cohrang.append(nois2coh[fc])
		nois3cohrang.append(nois3coh[fc])
		fc = fc + 14
	
	#Calculting the total effect of noise in the frequency domain
	#This stage keeps all the values to keep them at the same length as the next stage
	fillsnrlist = []
	fillsnrlist2 = []
	multsnrlist = []
	for mm in range(0,len(aveofnoiswin)):
		fillsnrlist.append(aveofnoiswin[mm]/sigwinpow[mm])
		fillsnrlist2.append(aveofnoiswin2[mm]/sigwinpow2[mm])
		multsnrlist.append((aveofnoiswin[mm]/sigwinpow[mm])*(aveofnoiswin2[mm]/sigwinpow2[mm]))
		
	#Keeping the range small for plotting
	finsnrs = []
	for snrc in range(0,len(multsnrs)):
		finsnrs.append(snrtracelist[snrc]+snrtracelist2[snrc]+multsnrs[snrc])
	
	
	#Saving the different nsrs for each stations
	s1n1.append(fillsnrlist) #These are not representative of the actual sums these would be n1/s1, n1n2/s1s2
	s2n2.append(fillsnrlist2)
	s1s2n1n2.append(multsnrlist)
	
	#Plotting the noise influence across frequencies for each station (s1n1+s2n2+s1s2n1n2)
	plt.figure(fignum)
	plt.xscale('log')
	plt.xlim((0.5,20))
	plt.title('Effect of noise in the frequency domain for station '+str(stats[j][0]))
	plt.ylabel('Noise contribution')
	plt.xlabel('Frequency (Hz)')
	plt.plot(frang, finsnrs,'b-')
	#bigplot, = plt.plot(frang, snrtracelist2,'r-')
	fignum = fignum+1
	
	#Storing the information
	allsnrslil.append(snrtracelist)
	allsnrsbig.append(snrtracelist2)
	sumsigcoh.append(sigcohrang)
	sumnois1coh.append(nois1cohrang)
	sumnois2coh.append(nois2cohrang)
	sumnois3coh.append(nois2cohrang)


#Also want to plot a sum of the nsrs across the frequency bands
#This may not be required, is left over from a previous version
finsnrslil = [sum(o) for o in zip(*allsnrslil)]
finsnrsbig = [sum(o) for o in zip(*allsnrsbig)]

#Averaging the real parts of the coherence between stations pairs for earthquake pair
#Using the five noise windows we have defined
finsigcoh = [sum(o)/len(o) for o in zip(*sumsigcoh)]
finnois1coh = [sum(o)/len(o) for o in zip(*sumnois1coh)]
finnois2coh = [sum(o)/len(o) for o in zip(*sumnois2coh)]
finnois3coh = [sum(o)/len(o) for o in zip(*sumnois3coh)]

#Plotting NSRs
#Again, it is unlikely this plot is of much use, but left in for my convenience
plt.figure(fignum)
plt.xscale('log')
plt.xlim((0.5,20))
plt.title('Summed Noise to Signal ratios in the frequency domain for all stations used')
plt.ylabel('Noise to Signal Ratio')
plt.xlabel('Frequency (Hz)')
aveplotlil, = plt.plot(frang, finsnrslil,'b-')
aveplotbig, = plt.plot(frang, finsnrsbig, 'r-')
plt.legend([aveplotlil,aveplotbig], ['M1.3 Earthquake', 'M3.7 Earthquake'])
fignum = fignum + 1

#Plotting phase coherence average
plt.figure(fignum)
finsigcohplot, = plt.plot(frang, finsigcoh, 'go')
#finnois1cohplot, = plt.plot(frang, finnois1coh, 'ro')
#finnois2cohplot, = plt.plot(frang, finnois2coh, 'co')
plt.xlim((0.5,20))
plt.title('Averaged phase coherence across all stations pairs for this specific earthquake pair')
plt.ylabel('Phase coherence')
plt.xlabel('Frequency (Hz)')
plt.legend([finsigcohplot], ['All Signal'])
#plt.legend([finsigcohplot,finnois1cohplot,finnois2cohplot,finnois3cohplot,finnois4cohplot], ['All Signal', 'All Presignal 1', 'All Presignal 2', 'All Postsignal 1', 'All Postsignal 2'])
fignum = fignum + 1
#So now we have several lists with unnormalised cross correlations within them, we want to then correlate between stations
#Setting a counter variable
n = 0

sumsig = []
sumnois1 = []
sumnois2 = []
sumnois3 = []
sumfnois = []
#Loop to perform the next level of correlations, and taking frequency spectra
#First loop goes through stations taken from list file earlier, second loop goes through the stations
#in that list after that station, so there are no repeats
for station in stats:
	for m in range(n+1,len(biglist)):
		#Assigning stations to be correlated, as well as pulling out the information we need to perform the next
		#correlation
		stat1 = stats[n][0]
		stat2 = stats[m][0]
		siggy1 = sigunnorm[n]
		siggy2 = sigunnorm[m]
		noisy1tr1 = nois1unnorm[n]
		noisy1tr2 = nois1unnorm[m]
		noisy2tr1 = nois2unnorm[n]
		noisy2tr2 = nois2unnorm[m]
		noisy3tr1 = nois3unnorm[n]
		noisy3tr2 = nois3unnorm[m]
		corrs1n1 = s1n1[n]
		corrs2n2 = s2n2[n] # These nsrs were taken in the frequency domain
		corrs3n3 = s1n1[m]
		corrs4n4 = s2n2[m]
		corrs1s2n1n2 = s1s2n1n2[n]
		corrs3s4n3n4 = s1s2n1n2[m]
		
		#Performing the unnormalised cross correlations between these different time windows
		fincorrsig = correlate(siggy1, siggy2, mode = 'full').tolist() #tolist() is required to allow these to work with prepwork
		fincorrnois1 = correlate(noisy1tr1, noisy1tr2, mode = 'full').tolist()
		fincorrnois2 = correlate(noisy2tr1, noisy2tr2, mode = 'full').tolist()
		fincorrnois3 = correlate(noisy3tr1, noisy3tr2, mode = 'full').tolist()

		#Fourier transforming and pulling out the frequency spectrum
		sigfft, sigconj, sigpow = prepwork(fincorrsig,1600)
		nois1fft, nois1conj, nois1pow = prepwork(fincorrnois1,1600)
		nois2fft, nois2conj, nois2pow = prepwork(fincorrnois2,1600)
		nois3fft, nois3conj, nois3pow = prepwork(fincorrnois3,1600)	
		
		#Calculating the total noise contribution
		stage1 = []
		stage2 = []
		stage3 = []
		fullnoisecont = []
		for nn in range(0,len(corrs1n1)):
			stage1.append(corrs1n1[nn]+corrs2n2[nn]+corrs3n3[nn]+corrs4n4[nn])
			stage2.append(corrs3s4n3n4[nn]+corrs1s2n1n2[nn]+(corrs1n1[nn]*corrs4n4[nn])+\
				(corrs2n2[nn]*corrs4n4[nn])+(corrs2n2[nn]*corrs3n3[nn])+(corrs1n1[nn]*corrs3n3[nn]))
			stage3.append((corrs3s4n3n4[nn]*corrs2n2[nn])+(corrs3s4n3n4[nn]*corrs1n1[nn])+\
				(corrs1s2n1n2[nn]*corrs3n3[nn])+(corrs1s2n1n2[nn]*corrs4n4[nn])+\
				(corrs1s2n1n2[nn]*corrs3s4n3n4[nn]))
			fullnoisecont.append(stage1[nn]+stage2[nn]+stage3[nn])
			
			
		#Taking only the frequency occurrence that we know we can sample (every 0.125Hz)
		sigrang = []
		nois1rang = []
		nois2rang = []
		nois3rang = []
		fnoisrang = []
		frang = []
		fc = 0
		num = len(fincorrsig)+1600
		timestep = 0.01
		freqs = np.fft.fftfreq(num, d = timestep)

		#Dividing by frequency as this transfers amplitudes to representative of displacement spectra
		#rather than velocity, as the traces were in originally
		while fc in range(0,3000):
			#fc = fc+3
			if fc != 0:
				frang.append(freqs[fc])
				sigrang.append(sigpow[fc]/freqs[fc])
				nois1rang.append(nois1pow[fc]/freqs[fc])
				nois2rang.append(nois2pow[fc]/freqs[fc])
				nois3rang.append(nois3pow[fc]/freqs[fc])
				fnoisrang.append(fullnoisecont[fc]*sigpow[fc])
				
			fc = fc + 14
		

		#Collecting values for averaging later
		#Plotting the noise and signal frequency spectrum for each station pair
		plt.figure(fignum+n)
		plt.subplot(len(biglist)-n-1,1,m-n)
		plt.yscale('log')
		plt.xscale('log')
		pairsigplot, = plt.plot(frang, sigrang, 'g-')
		#pairnois1plot, = plt.plot(frang, nois1rang, 'r-')
		#pairnois2plot, = plt.plot(frang, nois2rang, 'c-')
		#pairnois3plot, = plt.plot(frang, nois3rang, 'b-')
		fnoisranplot, = plt.plot(frang, fnoisrang, 'm-')
		plt.xlim((0.5,20))
		plt.title('Frequency spectra of signal and noise for '+stat1+' and '+stat2)
		plt.ylabel('Power')
		plt.xlabel('Frequency (Hz)')
		#plt.legend([pairsigplot,pairnois1plot,pairnois2plot,pairnois3plot, fnoisranplot], ['All Signal', 'Corr Presignal 1', 'Corr Presignal 2', 'Corr Presignal 3', 'NSR Noise signal'])
		#plt.legend([pairsigplot, fnoisranplot], ['All Signal',  'NSR Noise signal'])
		
		#Saving the values
		sumsig.append(sigrang)
		sumnois1.append(nois1rang)
		sumnois2.append(nois2rang)
		sumnois3.append(nois3rang)
		sumfnois.append(fnoisrang)
	
	#THIS IS VERY IMPORTANT, stations will not be looped through properly if this is not there
	n = n+1
	

#Summing all the spectra for the different time windows to maximise signal to noise
finsig = [sum(o) for o in zip(*sumsig)]
finnois1 = [sum(o) for o in zip(*sumnois1)]
finnois2 = [sum(o) for o in zip(*sumnois2)]
finnois3 = [sum(o) for o in zip(*sumnois3)]
finnsr = [sum(o) for o in zip(*sumfnois)]
fignum = fignum+n
#Plotting the sum of spectra
plt.figure(fignum)
plt.yscale('log')
plt.xscale('log')
sigsumplot, = plt.plot(frang, finsig, 'g-')
nois1sumplot, = plt.plot(frang, finnois1, 'r-')
nois2sumplot, = plt.plot(frang, finnois2, 'c-')
nois3sumplot, = plt.plot(frang, finnois3, 'b-')
fnoisplot, = plt.plot(frang, finnsr, 'm-')

plt.xlim((0.5,20))
plt.title('Frequency spectra of signal and noise for all station pairs summed')
plt.ylabel('Power')
plt.xlabel('Frequency (Hz)')
plt.legend([sigsumplot,nois1sumplot,nois2sumplot,nois3sumplot, fnoisplot], ['All Signal', 'Corr Presignal 1', 'Corr Presignal 2', 'Corr Presignal 3', 'Summed Noise results'])

plt.show()



