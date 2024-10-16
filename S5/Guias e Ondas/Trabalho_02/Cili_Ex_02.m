% Parâmetros do guia de onda cilíndrico
a = 0.025; % Raio (m)
c = 3e8;    % Velocidade da luz (m/s)
f = 12e9;   

% Primeiro zero da função de Bessel J_0 
x_01 = 2.405; % Zero da função de Bessel para TM01

% Frequência de corte
fc = (x_01 * c) / (2 * pi * a);
fprintf('Frequência de corte (TM01): %.2f GHz\n', fc / 1e9);

% Verificar se o modo é propagante
if f > fc
    fprintf('O modo TM01 é propagante.\n');
else
    fprintf('O modo TM01 não é propagante.\n');
end
%%
% Geração de uma malha de pontos cilíndricos
r = linspace(0, a, 100);
phi = linspace(0, 2*pi, 100);
[R, Phi] = meshgrid(r, phi);

% Campo elétrico E_z no modo TM01
Ez = besselj(0, x_01 * R / a);

% Conversão para coordenadas cartesianas
X = R .* cos(Phi);
Y = R .* sin(Phi);

% Visualização do campo elétrico
figure;
surf(X, Y, Ez);
title('Campo Elétrico (E_z) no Modo TM_{01}');
xlabel('x (m)');
ylabel('y (m)');
zlabel('E_z');
