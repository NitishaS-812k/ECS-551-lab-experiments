import numpy as np
from scipy import special
import matplotlib.pyplot as plt 

SNRdB = np.arange(-10,12,2)         #Signal to Noise ratio in dB
SNR = np.power(10, SNRdB/10)        #Convert SNR to normal scale
Arg_Erfc = np.sqrt(SNR)             #argument of erfc
BER_th = 0.5*(special.erfc(Arg_Erfc)) #theoretical BER
BER_sim = np.zeros(len(SNR))
max_run = 20
num_bit = 10000
Eb = 1
for i in range(len(SNR)):
    avgError = 0
    for _ in range(max_run):
        error = 0
        data = np.random.randint(2, size = num_bit)    #generate binary data
        X = 2*data - 1                                 # binary bpsk modulation
        mu = 0
        No = Eb/SNR[i]
        sigma = np.sqrt(No/2)
        N = sigma*np.random.randn(num_bit) + mu         #generate awgn channel
        Y = X + N                                       #received signal
        for k in range(num_bit):
            if (Y[k] > 0 and data[k] == 0) or (Y[k] < 0 and data[k] == 1):
                error += 1
        error = error/num_bit
        avgError = avgError + error
    BER_sim[i] = avgError/max_run

#plot BER
plt.plot(SNRdB, BER_th, 'r')
plt.plot(SNRdB, BER_sim, 'ko')
plt.legend(['Analytical', 'Simulation'])
plt.yscale('log')
plt.grid(True)
plt.ylabel('BER')
plt.xlabel('SNR(dB)')
plt.show()