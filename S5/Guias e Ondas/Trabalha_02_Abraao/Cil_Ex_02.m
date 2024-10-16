% Projeto Simples de Cálculo de Frequências de Corte em Guias de Onda Cilíndricos
clear; clc; close all;

% Parâmetros do guia de onda
c = 3e8; % Velocidade da luz no vácuo (m/s)

% Definição de raios para o guia de onda
raios = linspace(1e-2, 10e-2, 100); % Variando de 1 cm a 10 cm

% Zeros das funções de Bessel para os modos considerados
x_TE01 = 2.405; % Primeiro zero de J_0 para TE01
x_TE11 = 1.841; % Primeiro zero de J_1 para TE11
x_TM01 = 3.832; % Primeiro zero de J'_0 para TM01

% Inicialização dos vetores de frequências de corte
fc_TE01 = zeros(size(raios));
fc_TE11 = zeros(size(raios));
fc_TM01 = zeros(size(raios));

% Cálculo das frequências de corte para cada raio
for i = 1:length(raios)
    a = raios(i); % Raio do guia
    fc_TE01(i) = (x_TE01 * c) / (2 * pi * a); % Frequência de corte TE01
    fc_TE11(i) = (x_TE11 * c) / (2 * pi * a); % Frequência de corte TE11
    fc_TM01(i) = (x_TM01 * c) / (2 * pi * a); % Frequência de corte TM01
end

% Plot das frequências de corte
figure;
plot(raios, fc_TE01 * 1e-9, 'r', 'LineWidth', 2); % TE01
hold on;
plot(raios, fc_TE11 * 1e-9, 'b', 'LineWidth', 2); % TE11
plot(raios, fc_TM01 * 1e-9, 'g', 'LineWidth', 2); % TM01
hold off;

% Configurações do gráfico
xlabel('Raio do Guia (m)');
ylabel('Frequência de Corte (GHz)');
title('Frequências de Corte para Modos em Guias de Onda Cilíndricos');
legend('TE_{01}', 'TE_{11}', 'TM_{01}');
grid on;
