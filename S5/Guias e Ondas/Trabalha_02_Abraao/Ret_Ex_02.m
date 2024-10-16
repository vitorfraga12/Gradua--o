% Projeto de Simulação de Guias de Onda Retangulares em MATLAB
clear; clc; close all;

% Parâmetros do guia de onda retangular (Seção 1)
a1 = 5e-2; % Largura do guia (m)
b1 = 2e-2; % Altura do guia (m)

% Parâmetros do guia de onda retangular (Seção 2)
a2 = 4e-2; % Nova largura do guia (m)
b2 = 1.5e-2; % Nova altura do guia (m)

% Parâmetros do sinal
f = 10e9; % Frequência de operação (Hz)
c = 3e8; % Velocidade da luz no vácuo (m/s)
lambda = c / f; % Comprimento de onda (m)
omega = 2 * pi * f; % Frequência angular (rad/s)
k = 2 * pi / lambda; % Número de onda no vácuo (rad/m)

% Número de modos para cálculo
m = 1; % Modo na direção x
n = 0; % Modo na direção y

% Frequências de corte para as duas seções
fc1 = (c/2) * sqrt((m/a1)^2 + (n/b1)^2); % Frequência de corte (Seção 1)
fc2 = (c/2) * sqrt((m/a2)^2 + (n/b2)^2); % Frequência de corte (Seção 2)

% Verificação de propagação na Seção 1 e Seção 2
if f < fc1 || f < fc2
    error('A frequência de operação deve ser maior que a frequência de corte em ambas as seções.');
end

% Constantes de propagação
beta1 = sqrt(k^2 - (m*pi/a1)^2); % Constante de propagação Seção 1
beta2 = sqrt(k^2 - (m*pi/a2)^2); % Constante de propagação Seção 2

% Coeficiente de reflexão e transmissão na junção
Z1 = (k / beta1) * 377; % Impedância da Seção 1
Z2 = (k / beta2) * 377; % Impedância da Seção 2
Gamma = (Z2 - Z1) / (Z2 + Z1); % Coeficiente de reflexão
T = 1 + Gamma; % Coeficiente de transmissão

% Display dos resultados
fprintf('Frequência de corte na Seção 1: %.2f GHz\n', fc1*1e-9);
fprintf('Frequência de corte na Seção 2: %.2f GHz\n', fc2*1e-9);
fprintf('Constante de propagação na Seção 1: %.2f rad/m\n', beta1);
fprintf('Constante de propagação na Seção 2: %.2f rad/m\n', beta2);
fprintf('Coeficiente de Reflexão: %.2f\n', abs(Gamma));
fprintf('Coeficiente de Transmissão: %.2f\n\n', abs(T));

% Malha de pontos para visualização
x = linspace(0, 1, 200); % Direção de propagação
y = linspace(0, b1, 100); % Direção transversal
[X, Y] = meshgrid(x, y);

% Campo elétrico no guia (modo TE10)
E0 = 1; % Amplitude do campo elétrico
Ey1 = E0 * sin(m * pi * Y / b1) .* exp(-1j * beta1 * X); % Campo elétrico na Seção 1
Ey2 = E0 * sin(m * pi * Y / b2) .* exp(-1j * beta2 * X); % Campo elétrico na Seção 2

% Junção: Considerando descontinuidade na seção x = 0.5
Ey = Ey1; % Inicializa com a primeira seção
Ey(:, X(1, :) > 0.5) = Ey2(:, X(1, :) > 0.5); % Aplica descontinuidade

% Visualização do campo elétrico
figure;
imagesc(x, y, real(Ey));
colorbar;
title('Distribuição do Campo Elétrico (Re(E_y)) ao Longo do Guia de Onda');
xlabel('Direção de Propagação (x) [m]');
ylabel('Direção Transversal (y) [m]');

% Visualização do módulo do campo elétrico
figure;
imagesc(x, y, abs(Ey));
colorbar;
title('Módulo do Campo Elétrico |E_y| ao Longo do Guia de Onda');
xlabel('Direção de Propagação (x) [m]');
ylabel('Direção Transversal (y) [m]');
