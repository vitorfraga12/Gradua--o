from scipy import signal
import matplotlib.pyplot as plt

#ITEM 1
num1 = [1, 1]
den1 = [2, 3, 1]
sys1 = signal.TransferFunction(num1, den1)

# Calculando e plotando a resposta ao impulso
t1, y1 = signal.step(sys1)
t2, y2 = signal.impulse(sys1) 

plt.figure()
plt.plot(t1, y1)
plt.plot(t2, y2)
plt.legend(['Resposta ao Degrau', 'Resposta ao Impulso'])
plt.title('Resposta ao Degrau e ao Impulso para α = 1')
plt.xlabel('Tempo')
plt.ylabel('Amplitude')
plt.grid(True)
plt.show()

#ITEM 2.
alphas = [-1, -5, 0, 5, 7]
num = []   
den = [2, 3, 1]
# Criando a função de transferência
sys = []
for i in range(len(alphas)):
    num = [alphas[i], 1]
    sys.append(signal.TransferFunction(num, den))

t, y = [[] for _ in range(len(alphas))], [[] for _ in range(len(alphas))]


for i in range(len(alphas)):
    t[i], y[i] = signal.step(sys[i])    

plt.figure()
for i in range(len(alphas)):
    plt.plot(t[i], y[i])
plt.legend(['α = -1', 'α = -5', 'α = 0', 'α = 5', 'α = 7'])
plt.title('Resposta ao Degrau para diferentes valores de α')
plt.xlabel('Tempo')
plt.ylabel('Amplitude')
plt.grid(True)
plt.show()