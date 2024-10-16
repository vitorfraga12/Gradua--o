% Parâmetros
L = 0.05; % Comprimento da linha de carga (5 cm)
grid_size = 0.10; % Tamanho da grade (10 cm)
rho_cte = 1; % Constante de carga (ajuste conforme desejado)

% Constante de permitividade do vácuo
epsilon_0 = 8.854e-12; % F/m (farads por metro)

% Definir as coordenadas da grade tridimensional
x = linspace(-grid_size/2, grid_size/2, 100);
y = linspace(-grid_size/2, grid_size/2, 100);
z = linspace(-grid_size/2, grid_size/2, 100);

% Criar uma grade 3D
[X, Y, Z] = meshgrid(x, y, z);

% Calcular o campo elétrico (E) e o potencial elétrico (V)
E_x = zeros(size(X));
E_y = zeros(size(Y));
E_z = zeros(size(Z));
V = zeros(size(X));

for i = 1:numel(x)
    for j = 1:numel(y)
        for k = 1:numel(z)
            r = sqrt(X(i, j, k)^2 + Y(i, j, k)^2 + Z(i, j, k)^2);
            E_x(i, j, k) = (rho_cte * X(i, j, k)) / (2 * pi * epsilon_0 * r^3);
            E_y(i, j, k) = (rho_cte * Y(i, j, k)) / (2 * pi * epsilon_0 * r^3);
            E_z(i, j, k) = (rho_cte * Z(i, j, k)) / (2 * pi * epsilon_0 * r^3);
            V(i, j, k) = (rho_cte / (4 * pi * epsilon_0)) * (1 / r);
        end
    end
end

% Plote o campo elétrico e o potencial elétrico nos planos x = 0 e z = 0
figure;

subplot(2, 2, 1);
quiver(z, y, squeeze(E_z(:, 1, :)), squeeze(E_y(:, 1, :)));
title('Campo Elétrico em x = 0 (Plano y-z)');
xlabel('z (metros)');
ylabel('y (metros)');
axis equal;

subplot(2, 2, 2);
contourf(z, y, squeeze(V(:, 1, :)));
colorbar;
title('Potencial Elétrico em x = 0 (Plano y-z)');
xlabel('z (metros)');
ylabel('y (metros)');

subplot(2, 2, 3);
quiver(x, y, squeeze(E_x(:,:,1)), squeeze(E_y(:,:,1)));
title('Campo Elétrico em z = 0 (Plano x-y)');
xlabel('x (metros)');
ylabel('y (metros)');
axis equal;

subplot(2, 2, 4);
contourf(x, y, squeeze(V(:,:,1)));
colorbar;
title('Potencial Elétrico em z = 0 (Plano x-y)');
xlabel('x (metros)');
ylabel('y (metros)');

sgtitle('Campo e Potencial Elétrico (ρ = 1)');

%% 
% Parâmetros
L = 0.05; % Comprimento da linha de carga (5 cm)
grid_size = 0.10; % Tamanho da grade (10 cm)

% Constante de permitividade do vácuo
epsilon_0 = 8.854e-12; % F/m (farads por metro)

% Definir as coordenadas da grade tridimensional
x = linspace(-grid_size/2, grid_size/2, 100);
y = linspace(-grid_size/2, grid_size/2, 100);
z = linspace(-grid_size/2, grid_size/2, 100);

% Criar uma grade 3D
[X, Y, Z] = meshgrid(x, y, z);

% Calcular o campo elétrico (E) e o potencial elétrico (V) para rho = z
E_x = zeros(size(X));
E_y = zeros(size(Y));
E_z = zeros(size(Z));
V = zeros(size(X));

for i = 1:numel(x)
    for j = 1:numel(y)
        for k = 1:numel(z)
            r = sqrt(X(i, j, k)^2 + Y(i, j, k)^2 + Z(i, j, k)^2);
            rho = abs(Z(i, j, k)); % Distribuição de carga varia com |z|
            
            % Campo elétrico
            E_x(i, j, k) = (rho * X(i, j, k)) / (4 * pi * epsilon_0 * r^3);
            E_y(i, j, k) = (rho * Y(i, j, k)) / (4 * pi * epsilon_0 * r^3);
            E_z(i, j, k) = (rho * Z(i, j, k)) / (4 * pi * epsilon_0 * r^3);
            
            % Potencial elétrico
            V(i, j, k) = (rho / (4 * pi * epsilon_0)) * (1 / r);
        end
    end
end

% Plote o campo elétrico e o potencial elétrico nos planos x = 0 e z = 0
figure;

subplot(2, 2, 1);
quiver(z, y, squeeze(E_z(:, 1, :)), squeeze(E_y(:, 1, :)));
title('Campo Elétrico em x = 0 (Plano y-z)');
xlabel('z (metros)');
ylabel('y (metros)');
axis equal;

subplot(2, 2, 2);
contourf(z, y, squeeze(V(:, 1, :)));
colorbar;
title('Potencial Elétrico em x = 0 (Plano y-z)');
xlabel('z (metros)');
ylabel('y (metros)');

subplot(2, 2, 3);
quiver(x, y, squeeze(E_x(:,:,1)), squeeze(E_y(:,:,1)));
title('Campo Elétrico em z = 0 (Plano x-y)');
xlabel('x (metros)');
ylabel('y (metros)');
axis equal;

subplot(2, 2, 4);
contourf(x, y, squeeze(V(:,:,1)));
colorbar;
title('Potencial Elétrico em z = 0 (Plano x-y)');
xlabel('x (metros)');
ylabel('y (metros)');

sgtitle('Campo e Potencial Elétrico (ρ = z)');