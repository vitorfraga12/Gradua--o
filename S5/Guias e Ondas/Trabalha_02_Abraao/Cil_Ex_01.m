% Parâmetros do guia de onda cilíndrico
a = 5e-2; % Raio do guia de onda (m)
f = 15e9; % Frequência de operação (Hz)
c = 3e8; % Velocidade da luz no vácuo (m/s)
mu_0 = 4*pi*1e-7; % Permeabilidade do vácuo (H/m)
epsilon_0 = 8.854e-12; % Permissividade do vácuo (F/m)
eta = sqrt(mu_0/epsilon_0); % Impedância intrínseca do vácuo

% Cálculo da frequência de corte para o modo TE11
m = 1; % Ordem da função de Bessel (modo TE11)
n = 1; % Número de zeros da função de Bessel
x_mn = 1.841; % Primeiro zero da função de Bessel J_1
fc = (x_mn * c) / (2 * pi * a); % Frequência de corte (Hz)

% Verifica se a frequência de operação é maior que a frequência de corte
if f < fc
    error('A frequência de operação deve ser maior que a frequência de corte.');
end

% Cálculo do comprimento de onda no vácuo e número de onda
lambda = c / f; % Comprimento de onda no vácuo (m)
k = 2 * pi / lambda; % Número de onda no vácuo

% Constante de propagação no guia de onda
kc = x_mn / a; % Número de onda de corte
beta = sqrt(k^2 - kc^2); % Constante de propagação no guia

% Malha de pontos para visualização
rho = linspace(0, a, 100); % Direção radial
phi = linspace(0, 2*pi, 100); % Direção angular
[RHO, PHI] = meshgrid(rho, phi); % Malha em coordenadas cilíndricas

% Cálculo dos campos para o modo TE11
H0 = 1; % Amplitude do campo magnético
Hz = H0 * besselj(m, kc * RHO) .* cos(m * PHI); % Componente longitudinal do campo magnético
omega = 2 * pi * f; % Frequência angular
E_phi = -(H0 * beta / (omega * epsilon_0)) * besselj_prime(m, kc * RHO) .* sin(m * PHI); % Componente do campo elétrico transversal

% Visualização do campo magnético longitudinal Hz
figure;
polarplot3d(abs(Hz));
title('Distribuição do Campo Magnético Longitudinal |H_z| (Modo TE_{11})');
xlabel('Direção x (m)');
ylabel('Direção y (m)');
zlabel('|H_z| (A/m)');
colorbar;

% Visualização do campo elétrico transversal E_phi
figure;
polarplot3d(abs(E_phi));
title('Distribuição do Campo Elétrico Transversal |E_\phi| (Modo TE_{11})');
xlabel('Direção x (m)');
ylabel('Direção y (m)');
zlabel('|E_\phi| (V/m)');
colorbar;

% Função auxiliar para derivada de Bessel de primeira espécie
function J_prime = besselj_prime(m, x)
    J_prime = (besselj(m-1, x) - besselj(m+1, x)) / 2;
end

% Função para plotagem em coordenadas polares 3D
function polarplot3d(Z)
    [nRho, nPhi] = size(Z);
    rho = linspace(0, 1, nRho); % Normalização do raio
    phi = linspace(0, 2*pi, nPhi); % Ângulo
    [RHO, PHI] = meshgrid(rho, phi);
    X = RHO .* cos(PHI);
    Y = RHO .* sin(PHI);
    surf(X, Y, Z); % Sem parâmetros adicionais
    shading interp;
    axis equal tight;
    xlabel('x (m)');
    ylabel('y (m)');
    zlabel('|Z|');
end
