      SUBROUTINE phasecorr(master_trace, wavelet, npts, pcc)
 
      !------------------------------------------------------------------------
      ! This subroutine is called by a python script, which provides it with
      ! two time series containing the instananeous phase of two seismograms
      ! at each sample of the seismogram. It then performs the phase cross
      ! correlation method (equation 5 in Schimmel 1999), and returns the phase
      ! cross correlation between the two traces to the python script.
      ! 
      ! George Taylor, February 2015, University of Leeds
      !------------------------------------------------------------------------
 
      ! Variable declarations
      IMPLICIT NONE
      INTEGER lag_time
      INTEGER window_length
      INTEGER t_sample
      INTEGER samp_count
      INTEGER, INTENT(IN) :: npts
      REAL, INTENT(IN), DIMENSION(0:(npts-1)) :: master_trace
      REAL, INTENT(IN), DIMENSION(0:(npts-1)) :: wavelet
      REAL, DIMENSION (0:(npts-1)) :: corr_values
      REAL, DIMENSION ((-npts+1):(npts-1)), INTENT(OUT) :: pcc

      ! lag_time = t in PCC equation, represents the cross correlation lag time
      !
      ! window_length = N in the PCC equation, represents the no. of samples that
      !                 are being multiplied and summed to calculate the 
      !                 correlation value at a given lag time
      !
      ! t_sample = tau in PCC equation, time sample index in the wavelet and
      !            master trace that is being multiplied and appended to 
      !            corr_values
      !
      ! samp_count = Integer which keeps track of the entries in to the corr_value
      !              array, ensuring they go in to the right index
      !
      ! npts = Input argument that defines how many samples are contained in the
      !        master and wavelet traces
      !
      ! master_trace = Input argument containing the instantaneous phase at each
      !                sample for the seismogram that has been designated as the 
      !                'master'.
      ! 
      ! wavelet = Input argument containing the instantaneous phase for the 
      !           seismogram that has been designated as the 'secondary'
      !
      ! corr_values = Array containing the correlation values within the window,
      !               priot to summation and normalisation 
      !
      ! pcc = Output argument containing the full phase cross correlation
      !       function, returned to the calling script
      
      ! Initialise lag time to be the final sample of the master trace
      lag_time = npts - 1
      ! Begin the loop over all required lag times for the cross correlation
      DO
          ! Break from the loop once all possible lag time have been completed
          IF (lag_time .LT.(-npts+1)) EXIT
          ! Calculate the window_length for the calculation of the correlation 
          ! value
          !WRITE(*,*) 'Calculating for lag time', lag_time
          window_length = npts - abs(lag_time)
          ! Make sure the start of the corr_value window is placed correctly
          IF (lag_time .LT. 0) THEN
              t_sample = abs(lag_time)
          ELSE
              t_sample = 0
          ENDIF
          ! Reset the corr_values to zero to remove the effect of the previous 
          ! window
          corr_values = 0
          samp_count = 0 
          ! Calculate the correlation value of each sample of the master and
          ! wavelet traces that are in the current window
          DO
              ! Exit the loop if all samples in the current window are 
              IF (samp_count .GE. window_length) EXIT
              corr_values(samp_count) = (abs(cos((master_trace(lag_time+t_sample)&
                                      - wavelet(t_sample)) / 2.0)) - &
                                      abs(sin((master_trace(lag_time+t_sample) &
                                      - wavelet(t_sample)) /2.0)))
              t_sample = t_sample + 1
              samp_count = samp_count + 1
          END DO
          ! Calculate the final correlation value for the current lag time
          pcc(lag_time) = (1.0 / (window_length)) * sum(corr_values)
          lag_time = lag_time - 1
      ENDDO

      END SUBROUTINE phasecorr

         



      
      
