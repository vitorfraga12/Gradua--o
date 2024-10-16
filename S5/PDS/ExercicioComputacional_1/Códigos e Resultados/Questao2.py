import numpy as np
import matplotlib.pyplot as plt
import sounddevice as sd
import librosa
def dif_tras (x):
    n = len(x)
    y = np.zeros(n)
    for i in range(n):
        if i == 0:
            y[i] = x[i]
        else:
            y[i] = x[i] - x[i-1]
    return y

x, Fs = librosa.load('fala_sino.wav', sr=None)
y = dif_tras(x)

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
