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
