import numpy as np
import lal
import lalsimulation
from pycbc.waveform import get_td_waveform
from pycbc.filter import matched_filter
from pycbc.filter import make_frequency_series
from pycbc.types import TimeSeries, FrequencySeries
from scipy.signal import butter, filtfilt
from scipy.signal.windows import tukey

########################################################################
# FUNCTIONS TO BE SUPPLIED TO THE STUDENTS TO USE

def make_template(m1,m2,fs,T,inv_psd,s1=0.,s2=0.,f_low=20.,d=None,phic=None,tc=None):
    """
    This function generates a GW signal from a compact binary coalescence
    You must supply it with:
    m1 - the mass of the primary object (m1>m2) (units if solar mass)
    m2 - the mass of the secondary object (units of solar mass)
    fs - the sampling frequency (units of Hz)
    inv_psd - the inverse of the power spectral density (units of Hz)
    d [optional]- the distance to the source (units of megaparsecs)
    phic [optional] - the phase at coalescence (units of radians)
    tc [optional] - the time of coalescence from start (ubits of seconds)
    """
    N = int(T*fs)           # the total number of time samples
    dt = 1 / fs             # the sampling time (sec)
    df = 1.0/T              # the frequency resolution
    f_low = 20.0            # lowest frequency of waveform (Hz)
    if d is None:           # set the distance to 1 Mpc if not specified
        d = 1.0

    # make waveform
    n = N + 1
    while n>N/2:         # keep trying until the waveform safely fits inside the observation
        hp, _ = lalsimulation.SimInspiralChooseTDWaveform(
                    m1 * lal.MSUN_SI, m2 * lal.MSUN_SI,   # masses we comvert to SI units
                    0, 0, s1, 0, 0, s2,                     # spins we set to zero
                    d*1e6*lal.PC_SI,                      # the distance we convert to SI units
                    0, 0,            # inclination and phase
                    0, 0, 0,         # longAscNodes, eccentricity, and meanPerAno
                    dt,              # the time resolution
                    f_low,f_low,     # the min and reference frequencies
                    lal.CreateDict(),           # an empty placeholder
                    lalsimulation.IMRPhenomD)   # the waveform approximant 
        n = hp.data.data.size     # get the length of the generated waveform
        f_low += 1    # increment the low frequency start to reduce the signal length
        if f_low>100.:
            return False, False

    # pad the data and place the template in the middle 
    st = np.zeros(N)             # make an empty vector
    st[:n] = hp.data.data        # fill in the start of it with the signal
    idx = np.argmax(np.abs(st))  # find the time index with the maximum amplitude
    st = np.roll(st,-idx+N//2)   # roll the data to put the max amplitude in the centre

    # whiten the template using the given inverse PSD
    st_w = whiten(st,inv_psd,fs)

    # highpass the data - cut out all noise and signal power below 20 Hz 
    st_wb = highpass(st_w,20,fs)

    # put the template back and then shift the template to the requested time and phase
    f = df*np.arange(N//2 + 1)         # define frequency vector
    phase_t, phase_p = 0.0, 0.0        # initialise phase corrections
    if tc is not None:                 # if we set the time of coalescence
        st_wb = np.roll(st_wb,-N//2)   # undo the roll to the centre
        phase_t = 2.0j*np.pi*f*tc      # define phase correction for time shift
    if phic is not None:               # if we set the phase at coalescence
        phase_p = 1.0j*phic            # define phase shift

    # Fourier transform and apply the phase correction
    st_wbf = np.fft.rfft(st_wb)*np.exp(-phase_t - phase_p)
    st_wb = np.fft.irfft(st_wbf)    # inverse transform back to the time domain

    # make a time vector
    t = np.arange(int(fs*T))/fs

    # return template
    return t, st_wb


def make_template_pycbc(m1,m2,fs,T,inv_psd,s1=0.,s2=0.,f_low=20.,d=None,phic=None,tc=None):
    """
    This function generates a GW signal from a compact binary coalescence
    You must supply it with:
    m1 - the mass of the primary object (m1>m2) (units if solar mass)
    m2 - the mass of the secondary object (units of solar mass)
    fs - the sampling frequency (units of Hz)
    inv_psd - the inverse of the power spectral density (units of Hz)
    d [optional]- the distance to the source (units of megaparsecs)
    phic [optional] - the phase at coalescence (units of radians)
    tc [optional] - the time of coalescence from start (ubits of seconds)
    """
    N = int(T*fs)           # the total number of time samples
    dt = 1 / fs             # the sampling time (sec)
    df = 1.0/T              # the frequency resolution
    #f_low = 20.0            # lowest frequency of waveform (Hz)
    if d is None:           # set the distance to 1 Mpc if not specified
        d = 1.0

    # make waveform
    if True:
        hp, _ = get_td_waveform(approximant="IMRPhenomD",
                                mass1=m1, mass2=m2,
                                spin1z=s1, spin2z=s2,
                                delta_t=dt, f_lower=f_low)

        hp.resize(N)#int(2*(N-1)))
        st = np.array(hp.cyclic_time_shift(hp.start_time))

    # whiten the template using the given inverse PSD
    st_w = whiten(st,inv_psd,fs)

    # highpass the data - cut out all noise and signal power below 20 Hz 
    st_wb = highpass(st_w,20.,fs)

    # put the template back and then shift the template to the requested time and phase
    f = df*np.arange(N//2 + 1)         # define frequency vector
    phase_t, phase_p = 0.0, 0.0        # initialise phase corrections
    if tc is not None:                 # if we set the time of coalescence
        st_wb = np.roll(st_wb,-N//2)   # undo the roll to the centre
        phase_t = 2.0j*np.pi*f*tc      # define phase correction for time shift
    if phic is not None:               # if we set the phase at coalescence
        phase_p = 1.0j*phic            # define phase shift

    # Fourier transform and apply the phase correction
    st_wbf = np.fft.rfft(st_wb)*np.exp(-phase_t - phase_p) 
    st_wb = np.fft.irfft(st_wbf)    # inverse transform back to the time domain

    # make a time vector
    t = np.arange(int(fs*T))/fs

    # return template
    return t, st_wb

# function to whiten data
def whiten(x, inv_psd, fs):
    """
    This function whitened a timeseries by converting to the frequency domain
    and dividing through by the sqrt of the power spectral density.
    You must supply it with:
    x - input timeseries
    inv_psd - the inverse power spetrcal density (units of Hz)
    fs - the sampling frequency (units of Hz)  
    """

    Nt = x.size                     # the length of the timeseries
    dt = 1.0/fs                          # the sampling time
    freqs = np.fft.rfftfreq(Nt, dt)      # the frequency vectore for this timeseries

    # whitening: transform to freq domain, divide by sqrt(psd), then transform back, 
    # taking care to get normalization right.
    dwindow = tukey(Nt, alpha=1./16)             # define a tukey window
    hf = np.fft.rfft(x*dwindow)                  # apply Fourier transform to windowed data
    norm = 1./np.sqrt(1./(dt*2))                 # define normalisation
    white_hf = hf * np.sqrt(inv_psd) * norm      # devide through by sqrt(PSD) and normalised
    white_ht = np.fft.irfft(white_hf, n=Nt)      # inverse Fourier transform back to time domain
    
    # return whitened timeseries
    return white_ht

def highpass(x,flow,fs):
    """
    This function applies a 4th order high pass butterworth filter to a
    timeseries.
    You must supply it with:
    x - the input timeseries
    flow - the filter high pass frequency cut-off (units of Hz)
    fs - the sampling frequency (units if Hz)
    """

    # We need to suppress the high frequency noise (no signal!) with some bandpassing:
    bb, ab = butter(4, flow*2./fs, btype='highpass')    # make the filter
    normalization = np.sqrt((fs/2.0-flow)/(fs/2.0))     # compute normalisation
    return filtfilt(bb, ab, x) / normalization          # return after applying the filter and normalisation

def get_snr(data,template,fs):
    """
    Computes the SNR timeseries when given the whitened and normalised strain
    together with the whitened and normalised template
    You must supply
    data - the gravitational wave data supplied to you (dimensionless)
    template - a noise-free waveform that you have generated (dimensionless)
    fs - the sampling frequwency (units of Hz)
    """

    N = data.size  # the data length
    dt = 1.0/fs    # the sampling time (sec)
    T = N*dt       # the observation length (sec)
    df = 1.0/T     # the frequency resolution (Hz)
    
    # construct another identical template pi/2 out of phase
    template_fft = np.fft.rfft(template)     # go into the frequency domain
    template_quad = np.fft.irfft(-1.j*template_fft)   # multiply by -i and inverse Fourier transform back to the time domain

    # define the complex template - normalisation is arbitrary
    complex_template = (template + template_quad*1.j)
    complex_template_fft = np.fft.fft(complex_template)

    # Take the Fourier Transform (FFT) of the data and normalise
    data_fft = np.fft.fft(data) * np.sqrt(dt/2.0)

    # -- Calculate the matched filter output in the time domain:
    # Multiply the Fourier Space template and data.
    # Taking the Inverse Fourier Transform (IFFT) of the filter output puts it back in the time domain,
    # so the result will be plotted as a function of time off-set between the template and the data:
    optimal = data_fft * complex_template_fft.conjugate()
    optimal_time = 2.0*np.fft.ifft(optimal)*fs

    # -- Normalize the matched filter output:
    # Normalize the matched filter output so that we expect an average  value of 1 at times of just noise.
    # Then, the peak of the matched filter output will tell us the signal-to-noise ratio (SNR) of the signal.
    sigmasq = 1*(complex_template_fft * complex_template_fft.conjugate()).sum() * df
    sigma = np.sqrt(np.abs(sigmasq))
    SNR_complex = optimal_time/sigma

    # the template was set up to have peak amplitude in the middle of the tiome
    # vector so we now shift the SNR vector so that the peak SNR corresponds to
    # the actual time of the event
    SNR_complex = np.roll(SNR_complex,-N//2)
    SNR = abs(SNR_complex)  # take the absolute value of the complex SNR

    # return the SNR timeseries
    return SNR


