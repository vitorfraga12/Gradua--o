import numpy as np
from sympy import inverse_laplace_transform, symbols, Matrix, eye, Function, simplify, fraction, exp, DiracDelta, lambdify, laplace_transform
import matplotlib.pyplot as plt
from scipy import signal

#ITEM 1
#Sabemos que para encontrar a função de transferência dessa questão é necessário utilizar da seguinte equação:
#G(s) = C(Is - A)^-1 * B + D
#Declarando as variáveis de acordo com a equação

s, t = symbols('s t')

A = Matrix([[0, 1], [-5, -2]])
B = Matrix([[0], [2]])
C = Matrix([[0, 1]])
Is = s*eye(2)

G_s = (C * (Is - A).inv() * B)
G_s = G_s[0]
print(f'\nA função de transferência do sistema é:\n {G_s}')

#ITEM 2
#Temos que G(s) = Y(s)/U(s), onde Y(s) é a saída e U(s) é a entrada, fazendo então a separação das variáveis:
g_t = inverse_laplace_transform(G_s, s, t)
print(f'\nTemos então que o modelo entrada-saída do sistema é:\n {g_t}')

#ITEM 3
#Para encontrar a resposta a u(t) =  e^(−3t)*δ(t)
#Sabemos que G(s) = 2*s/(s**2 + 2*s + 5), quando multiplicamos por U(s) = 1/s + 3, temos que Y(s) = 2*s/(s**2 + 2*s + 5) * 1/(s+3)
#Logo Y(s) = 2*s/(s**3 + 5*s**2 + 11*s + 15)

num = [2, 0]
den = [1, 5, 11, 15]
sys = signal.TransferFunction(num, den)
t, y = signal.impulse(sys)

plt.figure()
plt.plot(t, y, label='y(t)')
plt.xlabel('Tempo (s)')
plt.ylabel('y(t)')
plt.title('Saída y(t) em resposta a u(t) = e^(-3t)*δ(t)')
plt.legend()
plt.grid(True)
plt.show()
