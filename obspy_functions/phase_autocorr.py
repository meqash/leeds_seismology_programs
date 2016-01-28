#!/usr/bin/env python

"""
This module contains functions for performing the phase autocorrelation
of ambient seismic noise after Schimmel et al. (2011). It includes
functions for performing the autocorrelation on one days worth of data,
as well as functions for stacking the autocorrelations.
George Taylor, University of Leeds
January 2016

"""

def autocorr(trace):
    """This function takes an obspy trace object and performs a phase autocorrelation
    of the trace with itself. First, the hilbert transform is taken to obtain the 
    analytic signal and hence the instantaneous phase. This is then passed to a 
    fortran script where phase correlation is performed after Schimmel et al., 2011.

    This function relies on the shared object file phasecorr.so, which is the file
    containing the fortran subroutine for phase correlation.
    """

    import numpy as np

    from scipy.signal import hilbert
    from phasecorr import phasecorr
    # Hilbert transform to obtain the analytic signal
    htrans = hilbert(trace)
    # Perform phase autocorrelation with instantaneous phase
    pac = phasecorr(np.angle(htrans),np.angle(htrans),len(htrans))

    return pac

def stream_stack(st):
    """
    This function takes an obspy stream object containing an arbitrary number of traces 
    and linearly stacks all traces contained within the stream.
    """

    import numpy as np

    return (np.sum([tr.data for tr in st], axis=0))/len(st)

def inst_correct(tr,pre_filt,unit):
    """This function corrects a trace from the DANA network for instrument
       response as read from the correct RESP file. It requires the trace
       object to be corrected, a pre filter (as a Python tuple with 4 entries)
       and a string object stipulating the unit to be corrected to (either
       "DIS", "VEL" or "ACC")
    """
    import os
    
    curr_dir = os.getcwd() # Remember working directory

    resp_dir = "/nfs/a224/DANA/DATALESS/RESPS" # Directory containing RESPS
    os.chdir(resp_dir)

    respf = "RESP.YH."+tr.stats.station+".."+tr.stats.channel
    seedresp = {"filename": respf,
                 "units": unit
                }

    corr_tr = tr.simulate(paz_remove=None, pre_filt=pre_filt, seedresp=seedresp)
    os.chdir(curr_dir) # Change back to previous directory
    return corr_tr

def mktrace(data,station,channel,sampling_rate,npts):
    """This function creates a new obspy trace object given the array containing
       trace data and important meta data.

       It requires a numpy array containing the data to be turned into a trace,
       and string objects of the station name and channel, and integer values
       stipulating the sampling rate and total number of points in the data.

       An obspy Trace object is returned.
    """
    from obspy.core import Trace

    stats={"station": station,
           "channel": channel,
           "sampling_rate": sampling_rate,
           "npts": npts
          }
    return Trace(data=data,header=stats)

def gen_fname(name,net,chan,year,jul_day):
    """This function takes trace information including: station name, network,
    channel, year and julian day. This information is then converted into a
    useful filename for the trace.
    """
    return str(name+"."+net+"."+chan+"."+str(year)+"."+str(jul_day))
