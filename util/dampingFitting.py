# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
#    p    #     Version: 0.4
#    y    #     Date: 03/11/2023
#    F    #     Author: Martin Saravia
#    S    #     Description: Function to calculate damping coefficient of a decaying signal
#    I    #     Return: Damping parameters
# --------------------------------------------------------------------------- #
# Notes:
#   This functions calculate the decay rate of a signal. 
#
# Warnings:
#   
# Optimize:
#   
# --------------------------------------------------------------------------- #

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from scipy.optimize import curve_fit
from scipy.signal import butter
from scipy.signal import sosfilt

def dampingFit(time, signal, frequency, cutEnvelope, filter=False):

    # Calculate the sampling rate from the time vector
    samplingRate = int(1 / (time[1] - time[0]))
    print("The sampling rate of the signal is: ", samplingRate)
    
    # Center the signal
    signal -= np.mean(signal)
    if filter:
        # Filter the signal using a Butterworth filter
        sos = butter(10, filter, btype='bandpass', fs=samplingRate, output='sos')

        signal = sosfilt(sos, signal)

    # Estimated natural frequency of the signal (we need it to get the proportional damping ratio)
    omega = frequency * 6.28 

    cut0 = 0 # start point of the envelope we want to process (choose is to avoid zeros)
    cut1 = -1 # end point of the envelope

    # Hilbert transform of the signal
    Ht = hilbert(signal)

    # Extract the envelope by calculating the magnitude of the original signal
    envelope = np.abs(Ht)

    # Plot the original signal and its Hilbert transform
    plt.figure(figsize=(10, 6))
    plt.subplot(3, 1, 1)
    plt.plot(time, np.real(signal), label='Real Part')
    plt.plot(time, np.imag(signal), label='Imaginary Part')
    plt.title('Complex Exponential Signal')
    plt.legend()

    plt.subplot(3, 1, 2)
    plt.plot(time, np.real(Ht), label='Hilbert Transform (Real Part)')
    plt.plot(time, np.imag(Ht), label='Hilbert Transform (Imaginary Part)')
    plt.title('Hilbert Transform')
    plt.legend()



    # Exponential decay function
    def damped_exponential(t, A, decay_constant, C):
        return A * np.exp(-decay_constant * t) + C

    # Fit the amplitude envelope to a exponential decay function
    timeCut = cutEnvelope / samplingRate 
    fitTime = time[cut0:cut1] - timeCut # We MUST shift the time to fit
    fitEnvelope = envelope[cut0:cut1]
    p0 = (np.max(fitEnvelope), 30, np.min(fitEnvelope)) # Initial guess

    print("Fitting with initial guess: ", p0)
    params, _ = curve_fit(damped_exponential, fitTime, fitEnvelope, p0=p0)
    print("The polynomial fitted damping ratio is: ", params[1]/omega)
    print("The polynomial fitted parameters are: ", params)

    plt.subplot(3, 1, 3)
    plt.plot(fitTime, fitEnvelope, label='Hilbert Envelope')
    fit = damped_exponential(fitTime, params[0], params[1], params[2])
    plt.plot(fitTime, fit, label='Exponential fit')
    plt.title('Hilbert fitting')
    plt.legend()

    plt.tight_layout()
    plt.show()