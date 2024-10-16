import numpy as np
import matplotlib.pyplot as plt
import sounddevice as sd
import librosa
def med_mov (x, M1, M2):
    n = len(x)
    y = np.zeros(n)
    for i in range(n):
        sum_x = 0
        for k in range(-M1, M2+1):
            if i+k >= 0 and i+k < n:
                sum_x += x[i+k]
        y[i] = sum_x/(M1+M2+1)
    return y

x, Fs = librosa.load('fala_sino.wav', sr=None)
M1, M2 = 5, 10
y = med_mov(x, M1, M2)

#Análise no Tempo
plt.figure()
plt.plot(x)
plt.plot(y)
plt.legend(['Aúdio Original', 'Aúdio Filtrado'])
plt.xlabel('Tempo')
plt.show()

#Análise na Frequência
X = np.fft.fft(x)
Y = np.fft.fft(y)

# Calcula as frequências correspondentes
freq = np.fft.fftfreq(len(x), 1/Fs)

X = X[:len(X)//2]
Y = Y[:len(Y)//2]
freq = freq[:len(freq)//2]

# Gráfico na frequência
plt.figure()
plt.plot(freq, np.abs(X))
plt.plot(freq, np.abs(Y))
plt.legend(['FFT Original', 'FFT Filtrado'])
plt.xlabel('Frequência (Hz)')
plt.show()

# Rodar o áudio
sd.play(x, Fs)
sd.wait()
sd.play(y, Fs)
sd.wait()

