% Definindo parâmetros
k = 8.99e9;  % Constante eletrostática em N.m^2/C^2
a = 1;      % Definindo o raio do disco
[X, Y] = meshgrid(linspace(-a, a, 100), linspace(-a, a, 100));  % Malha de pontos x e y

% Convertendo x e y para rho
Rho = sqrt(X.^2 + Y.^2);

% Evitar divisão por zero no centro
Rho(Rho < 1e-10) = 1e-10;

% Campo Elétrico
E = arrayfun(@(r) integral(@(rp) k*rp.^2./r.^2, 0, a), Rho);

% Potencial Elétrico (Uso da integral cumulativa trapezoidal)
V = -cumtrapz(Rho(:), E(:));
V = reshape(V, size(Rho));

% Plotando
figure;

subplot(1, 2, 1);
surf(X, Y, E, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
zlabel('Campo Elétrico (E)');
title('Campo Elétrico');

subplot(1, 2, 2);
surf(X, Y, V, 'EdgeColor', 'none');
colorbar;
xlabel('x');
ylabel('y');
zlabel('Potencial Elétrico (V)');
title('Potencial Elétrico');