% Parâmetros do guia cilíndrico
a = 0.025; % raio (m)
c = 3e8;    % velocidade da luz (m/s)
f = 10e9;   % Frequência

% Função de Bessel e zeros
n = 1; % ordem do modo TE
m = 1; % primeiro zero da função de Bessel J_n

% Primeiro zero da função de Bessel J_n
x_nm = 2.405; % primeiro zero de J_1 (modo TE11)

% Frequência de corte
fc = (x_nm * c) / (2 * pi * a);

% Constante de propagação
lambda_c = c / fc;
lambda = c / f;
beta = 2 * pi * sqrt(1/lambda^2 - 1/lambda_c^2);

fprintf('Frequência de corte: %.2f GHz\n', fc / 1e9);
fprintf('Constante de propagação β: %.2f rad/m\n', beta);

% Geração de uma malha de pontos em coordenadas cilíndricas
r = linspace(0, a, 100);
phi = linspace(0, 2*pi, 100);
[R, Phi] = meshgrid(r, phi);

% Função de Bessel e campo magnético H_z no modo TE11
Hz = besselj(n, x_nm * R / a) .* cos(n * Phi);

% Conversão para coordenadas cartesianas
X = R .* cos(Phi);
Y = R .* sin(Phi);

% Visualização do campo magnético
figure;
surf(X, Y, Hz);
title('Campo Magnético (H_z) no Modo TE_{11}');
xlabel('x (m)');
ylabel('y (m)');
zlabel('H_z');
