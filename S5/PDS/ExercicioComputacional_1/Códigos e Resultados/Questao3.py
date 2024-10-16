import numpy as np
import matplotlib.pyplot as plt
import sounddevice as sd
import librosa

def convolucao(x, h):
    N1 = len(x)
    N2 = len(h)
    N = N1 + N2 - 1
    y = np.zeros(N)
    for n in range(N):
        for k in range(N2):
            if n - k >= 0 and n - k < N1:
                y[n] += x[n-k] * h[k]
    return y

x, Fs = librosa.load('fala_sino.wav', sr=None)

n = len(x)
h = [1, 0.5, 0, -0.25]

y = convolucao(x, h)

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