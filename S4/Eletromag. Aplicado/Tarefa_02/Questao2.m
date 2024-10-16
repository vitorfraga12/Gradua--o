
% Parâmetros
a = 0.02; % Raio da espira em metros (2 cm)
I_constante = 1; % Corrente em ampères (constante)
grid_size = 0.10; % Tamanho da grade (10 cm)

% Constante permeabilidade do vácuo
mu_0 = 4 * pi * 1e-7; % H/m (henrys por metro)

% Definir as coordenadas da grade tridimensional
x = linspace(-grid_size/2, grid_size/2, 100);
y = linspace(-grid_size/2, grid_size/2, 100);
z = linspace(-grid_size/2, grid_size/2, 100);

% Criar uma grade 3D
[X, Y, Z] = meshgrid(x, y, z);

% Inicializar o campo magnético B nas direções x, y e z
B_x = zeros(size(X));
B_y = zeros(size(Y));
B_z = zeros(size(Z));

% Loop sobre os pontos da espira para I constante
for phi = linspace(0, 2 * pi, 100)
    x_spira = a * cos(phi);
    y_spira = a * sin(phi);
    z_spira = 0; % Plano z = 0
    
    % Vetor de posição do ponto da espira ao ponto onde queremos calcular B
    r = sqrt((X - x_spira).^2 + (Y - y_spira).^2 + (Z - z_spira).^2);
    
    % Componentes do campo magnético dBr, dBphi e dBz para I constante
    dBr = mu_0 * I_constante * a^2 * sin(phi) ./ (2 * pi * r.^3);
    dBphi = mu_0 * I_constante * a^2 * cos(phi) ./ (2 * pi * r.^3);
    dBz = mu_0 * I_constante * a^2 ./ (2 * pi * r.^3);
    
    % Somar as contribuições de cada ponto da espira ao campo magnético total
    B_x = B_x + dBr;
    B_y = B_y + dBphi;
    B_z = B_z + dBz;
end

% Plote o campo magnético no plano x = 0 (y-z) para I constante
figure;

subplot(1, 2, 1);
quiver(z, y, squeeze(B_z(:, 1, :)), squeeze(B_y(:, 1, :)));
title('Campo Magnético em x = 0 (Plano y-z) para I = 1');
xlabel('z (metros)');
ylabel('y (metros)');
axis equal;

% Plote o campo magnético no plano z = 0 (x-y) para I constante
subplot(1, 2, 2);
quiver(x, y, squeeze(B_x(:,:,1)), squeeze(B_y(:,:,1)));
title('Campo Magnético em z = 0 (Plano x-y) para I = 1');
xlabel('x (metros)');
ylabel('y (metros)');
axis equal;

sgtitle('Campo Magnético (I = 1)');

% Repita o cálculo para I = sen(φ) A
B_x_sen = zeros(size(X));
B_y_sen = zeros(size(Y));
B_z_sen = zeros(size(Z));

for phi = linspace(0, 2 * pi, 100)
    x_spira = a * cos(phi);
    y_spira = a * sin(phi);
    z_spira = 0; % Plano z = 0
    
    % Vetor de posição do ponto da espira ao ponto onde queremos calcular B
    r = sqrt((X - x_spira).^2 + (Y - y_spira).^2 + (Z - z_spira).^2);
    
    % Componentes do campo magnético dBr, dBphi e dBz para I = sen(φ)
    dBr = mu_0 * a^2 * sin(phi) .* sin(phi) ./ (2 * pi * r.^3);
    dBphi = mu_0 * a^2 * cos(phi) .* sin(phi) ./ (2 * pi * r.^3);
    dBz = mu_0 * a^2 * sin(phi) ./ (2 * pi * r.^3);
    
    % Somar as contribuições de cada ponto da espira ao campo magnético total para I = sen(φ)
    B_x_sen = B_x_sen + dBr;
    B_y_sen = B_y_sen + dBphi;
    B_z_sen = B_z_sen + dBz;
end

% Plote o campo magnético no plano x = 0 (y-z) para I = sen(φ)
figure;

subplot(1, 2, 1);
quiver(z, y, squeeze(B_z_sen(:, 1, :)), squeeze(B_y_sen(:, 1, :)));
title('Campo Magnético em x = 0 (Plano y-z) para I = sen(φ)');
xlabel('z (metros)');
ylabel('y (metros)');
axis equal;

% Plote o campo magnético no plano z = 0 (x-y) para I = sen(φ)
subplot(1, 2, 2);
quiver(x, y, squeeze(B_x_sen(:,:,1)), squeeze(B_y_sen(:,:,1)));
title('Campo Magnético em z = 0 (Plano x-y) para I = sen(φ)');
xlabel('x (metros)');
ylabel('y (metros)');
axis equal;

sgtitle('Campo Magnético (I = sen(φ))');
%% 
% Parâmetros
a = 0.02; % Raio da espira em metros (2 cm)
I_constante = 1; % Corrente em ampères (constante)
grid_size = 0.10; % Tamanho da grade (10 cm)

% Constante permeabilidade do vácuo
mu_0 = 4 * pi * 1e-7; % H/m (henrys por metro)

% Definir as coordenadas da grade tridimensional
x = linspace(-grid_size/2, grid_size/2, 100);
y = linspace(-grid_size/2, grid_size/2, 100);
z = linspace(-grid_size/2, grid_size/2, 100);

% Criar uma grade 3D
[X, Y, Z] = meshgrid(x, y, z);

% Inicializar o campo magnético B nas direções x, y e z
B_x = zeros(size(X));
B_y = zeros(size(Y));
B_z = zeros(size(Z));

% Loop sobre os pontos da espira para I constante
for phi = linspace(0, 2 * pi, 100)
    x_spira = a * cos(phi);
    y_spira = a * sin(phi);
    z_spira = 0; % Plano z = 0
    
    % Vetor de posição do ponto da espira ao ponto onde queremos calcular B
    r = sqrt((X - x_spira).^2 + (Y - y_spira).^2 + (Z - z_spira).^2);
    
    % Componentes do campo magnético dBr, dBphi e dBz para I constante
    dBr = mu_0 * I_constante * a^2 * sin(phi) ./ (2 * pi * r.^3);
    dBphi = mu_0 * I_constante * a^2 * cos(phi) ./ (2 * pi * r.^3);
    dBz = mu_0 * I_constante * a^2 ./ (2 * pi * r.^3);
    
    % Somar as contribuições de cada ponto da espira ao campo magnético total
    B_x = B_x + dBr;
    B_y = B_y + dBphi;
    B_z = B_z + dBz;
end

% Passo para espaçamento de vetores (aumente para reduzir o número de vetores)
step = 5;

% Plote o campo magnético no plano x = 0 (y-z) para I constante
figure;

subplot(1, 2, 1);
quiver(z(1:step:end, 1:step:end), y(1:step:end, 1:step:end), ...
    squeeze(B_z(1:step:end, 1:step:end, 1:step:end)), ...
    squeeze(B_y(1:step:end, 1:step:end, 1:step:end)));
title('Campo Magnético em x = 0 (Plano y-z) para I constante');
xlabel('z (metros)');
ylabel('y (metros)');
axis equal;

% Plote o campo magnético no plano z = 0 (x-y) para I constante
subplot(1, 2, 2);
quiver(x(1:step:end, 1:step:end), y(1:step:end, 1:step:end), ...
    squeeze(B_x(1:step:end, 1:step:end, 1)), ...
    squeeze(B_y(1:step:end, 1:step:end, 1)));
title('Campo Magnético em z = 0 (Plano x-y) para I constante');
xlabel('x (metros)');
ylabel('y (metros)');
axis equal;

sgtitle('Campo Magnético (I constante)');

% Repita o cálculo para I = sen(φ) A
B_x_sen = zeros(size(X));
B_y_sen = zeros(size(Y));
B_z_sen = zeros(size(Z));

for phi = linspace(0, 2 * pi, 100)
    x_spira = a * cos(phi);
    y_spira = a * sin(phi);
    z_spira = 0; % Plano z = 0
    
    % Vetor de posição do ponto da espira ao ponto onde queremos calcular B
    r = sqrt((X - x_spira).^2 + (Y - y_spira).^2 + (Z - z_spira).^2);
    
    % Componentes do campo magnético dBr, dBphi e dBz para I = sen(φ)
    dBr = mu_0 * a^2 * sin(phi) .* sin(phi) ./ (2 * pi * r.^3);
    dBphi = mu_0 * a^2 * cos(phi) .* sin(phi) ./ (2 * pi * r.^3);
    dBz = mu_0 * a^2 * sin(phi) ./ (2 * pi * r.^3);
    
    % Somar as contribuições de cada ponto da espira ao campo magnético total para I = sen(φ)
    B_x_sen = B_x_sen + dBr;
    B_y_sen = B_y_sen + dBphi;
    B_z_sen = B_z_sen + dBz;
end

% Plote o campo magnético no plano x = 0 (y-z) para I = sen(φ)
figure;

subplot(1, 2, 1);
quiver(z(1:step:end, 1:step:end), y(1:step:end, 1:step:end), ...
    squeeze(B_z_sen(1:step:end, 1:step:end, 1:step:end)), ...
    squeeze(B_y_sen(1:step:end, 1:step:end, 1:step:end)));
title('Campo Magnético em x = 0 (Plano y-z) para I = sen(φ)');
xlabel('z (metros)');
ylabel('y (metros)');
axis equal;

% Plote o campo magnético no plano z = 0 (x-y) para I = sen(φ)
subplot(1, 2, 2);
quiver(x(1:step:end, 1:step:end), y(1:step:end, 1:step:end), ...
    squeeze(B_x_sen(1:step:end, 1:step:end, 1)), ...
    squeeze(B_y_sen(1:step:end, 1:step:end, 1)));
title('Campo Magnético em z = 0 (Plano x-y) para I = sen(φ)');
xlabel('x (metros)');
ylabel('y (metros)');
axis equal;

sgtitle('Campo Magnético (I = sen(φ))');
