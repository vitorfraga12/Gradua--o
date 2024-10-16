% Parâmetros do guia de onda
l = 0.025; % largura (m)
h = 0.01; % altura (m)
c = 3e8;    % velocidade da luz (m/s)

% Definir a frequência de operação
f = 10e9;   %

% Cálculo da frequência de corte para modos TE_{10}
m = 1; % 
n = 0;

% Fórmula da frequência de corte
fc = (c / 2) * sqrt((m / l)^2 + (n / h)^2);

% Constante de propagação
lambda_c = c / fc;  % comprimento de onda de corte
lambda = c / f;     % comprimento de onda de operação
beta = 2 * pi * sqrt(1/lambda^2 - 1/lambda_c^2); % constante de propagação

fprintf('Frequência de corte: %.2f GHz\n', fc / 1e9);
fprintf('Constante de propagação β: %.2f rad/m\n', beta);

% Geração de uma malha para plotar os campos
x = linspace(0, l, 100);
y = linspace(0, h, 100);
[X, Y] = meshgrid(x, y);

% Equação do campo elétrico no modo TE10
Ez = sin(pi * X / l);

% Visualização do campo
figure;
surf(X, Y, Ez);
title('Campo Elétrico (E_z) no Modo TE_{10}');
xlabel('x (m)');
ylabel('y (m)');
zlabel('E_z');
