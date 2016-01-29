#specext appends zeroes, rotates, takes the fourier transform and extracts the power spectrum
#as well as the conjugate, it can also be adjusted to print out the phase by uncommenting out certain lines
#Two things must be input here, your time series (vals) and your length of time (tw - needs to be number of samples, not in seconds)
#that you want to rotate by. 
#Due to lengthening trace to 10000 by appending zeroes, frequency spacing will be every 0.01. Be aware you can change this number, 
#depending on how long the data series you expect to work with are
def specext(vals,tw):
	import numpy as np
	#This line appends zeroes onto the end of the trace, set at 10000 so the frequency information you get out will always be in the same
	#context, without this, you may be plotting the wrong frequency range with the wrong spectrum
	vals.extend([0]*(10000-len(vals)))
	
	#These two lines allow you to rotate the trace by a set time window which is specified in the function, this is useful for shifting
	#to the arrival time of an earthquake for example
	#A deque is a certain form of representing the data that allows you to do this
	valsmov = deque(vals)
	valsmov.rotate(-tw)
	
	#Following this, you take the fourier transform, followed by pulling out the power spectrum and the conjugate
	valsfft = np.fft.fft(valsmov)
	valspow = np.abs(valsfft)
	conjugate = np.conj(valsfft)
	#phas = np.angle(valsfft,deg=True)
	
	#Returning values for use
	return valsfft, conjugate, valspow #,phas



#sleptap trims a trace (in this case a copy of the original trorig), tapers it using a slepian taper (multi taper), and outputs the trace as a numpy array
#As well as outputting the array as a list (slepwin) for use in specext
#For this to function, tror must have been imported using obspy, your low time and hitime are the times you want to trim the trace down to
#slepper and slepwid define the percentage to be tapered of the window you've selected and the bandwidth to be used in the slepian taper
def sleptap(tror, lowtime, hitime, slepper, slepwid):
	#First pull out the start time for use in the function
	start = tror[0].stats.starttime
	
	#Copying the trace in order to keep the original version, it is then trimmed to the specified window
	sleptr = tror.copy().trim(start+(lowtime*0.01), start+(hitime*0.01))
	
	#Next we perform the taper using a slepian window, in this case both sides are tapered, but this can easily be changed by adjusting the 
	#side variable
	sleptr.taper(max_percentage=slepper,type='slepian',max_length=None,side='both',width=slepwid)
	
	#Pulling the trace out into a numpy array
	sleparr = sleptr[0]
	#Also outputting the trace in a list form
	slepwin = sleparr[:].tolist()
	
	#Returning values
	return sleparr, slepwin
