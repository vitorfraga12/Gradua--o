% MATLAB Code: Padrão de radiação para d = lambda/4

% Parâmetros iniciais
N = 4;                % Número de elementos no array
d = 0.25;             % Espaçamento normalizado (d/lambda)
theta = linspace(0, 2*pi, 1000); % Angulos de varredura [rad]
k = 2*pi;             % Constante de propagação

% Definição do fator de array
AF = zeros(size(theta));
for n = 1:N
    AF = AF + exp(1j*(n-1)*k*d*cos(theta));
end

% Normalização da magnitude do fator de array
AF_magnitude = abs(AF) / max(abs(AF));

% Gráfico polar
figure;
polarplot(theta, 20*log10(AF_magnitude), 'LineWidth', 1.5);
rlim([-40 0]); % Limite do eixo radial em dB
title('Padrão de Radiação para d = \lambda/4');

grid on;
set(gca, 'ThetaZeroLocation', 'top', 'ThetaDir', 'clockwise');

% Salvando o gráfico
saveas(gcf, 'AF_lambda4.png');
%%
% Intervalo para 0 <= theta <= 90 graus
theta = linspace(-pi/2, pi/2, 1000); % Primeiro quadrante (0 a 90 graus)
z_visible = exp(1j * (k * d * cos(theta) + beta)); % Região visível no 1º quadrante

% Gráfico ajustado
figure;
hold on;
plot(real(z_circle), imag(z_circle), 'k--', 'LineWidth', 1, 'DisplayName', 'Círculo Unitário Completo');
plot(real(z_visible), imag(z_visible), 'b', 'LineWidth', 2, 'DisplayName', 'Região Visível (0^\\circ \\leq \\theta \\leq 90^\\circ)');
scatter(real(z_visible(1)), imag(z_visible(1)), 50, 'r', 'filled', 'DisplayName', 'Início: \\theta = 0^\\circ');
scatter(real(z_visible(end)), imag(z_visible(end)), 50, 'g', 'filled', 'DisplayName', 'Fim: \\theta = 90^\\circ');

% Configurações do gráfico
grid on;
axis equal;
xlabel('Re(z) (Parte Real)');
ylabel('Im(z) (Parte Imaginária)');
title('Círculo Unitário e Região Visível no Primeiro Quadrante');
legend show;

% Salvando o gráfico
saveas(gcf, 'circulo_unitario_quadrante1.png');
%%
% Configuração de parâmetros
N = 4; % Número de elementos no array
d = 0.25; % Espaçamento normalizado (d/lambda)
theta_nulos = [30, 90, 150]; % Direções dos nulos em graus
k = 2*pi; % Constante de propagação
theta = linspace(0, 180, 1000); % Ângulos de varredura em graus

% Conversão de ângulos para radianos
psi_nulos = cosd(theta_nulos); % Projeção angular
z_nulos = exp(1j * pi * psi_nulos); % Cálculo das raízes (nulos)

% Fator de array usando o método de Schelkunoff
AF = ones(size(theta)); % Inicializa o fator de array
for zn = z_nulos
    AF = AF .* (exp(1j * pi * cosd(theta)) - zn); % Atualiza o fator de array
end

% Normalização do padrão de radiação
AF_magnitude = abs(AF) / max(abs(AF)); % Normaliza a magnitude do fator de array

% Gráfico polar do padrão de radiação
figure;
polarplot(deg2rad(theta), 20*log10(AF_magnitude), 'LineWidth', 1.5); % Gráfico polar
title('Padrão de Radiação com Nulos em 30^\circ, 90^\circ e 150^\circ'); % Título
rlim([-40 0]); % Limite radial em dB
grid on;

% Salvando o gráfico
saveas(gcf, 'padrao_radiacao_nulos_N4.png'); % Salva o gráfico como imagem



