from sympy import symbols, inverse_laplace_transform, exp

# Definindo as variáveis "t" no domínio do Tempo e "s" no domínio de Laplace
s, t = symbols('s t')

# Definindo as funções que eu quero encontrar a transformada inversa de Laplace
F1 = (4*s**2 + 15*s+ 19)/(s**2 + 2*s + 1 )
F2 = (s**2 + 2*s +1)/((s + 2)**3)
F3 = (2*s + 3)/(s**3 + 6*s**2 + 21*s + 26)
F4 = (1 + 2*exp(-s))/(s**2 + 3*s + 2)
# Encontrando as transformada inversas de Laplace
f1 = inverse_laplace_transform(F1, s, t)
f2 = inverse_laplace_transform(F2, s, t)
f3 = inverse_laplace_transform(F3, s, t)
f4 = inverse_laplace_transform(F4, s, t)

# Imprimir o resultado
print(f'Temos que as transformadas inversas de Laplace são:\n f1(t) = {f1} \n f2(t) = {f2} \n f3(t) = {f3} \n f4(t) = {f4}')


