%% %variação da capacitância com x0 
% Intervalo de deslocamento x0 de 0 a 5 cm
x0_min = 0; % Deslocamento mínimo em metros
x0_max = 0.05; % Deslocamento máximo em metros
num_pontos = 50; % Número de pontos no intervalo

% Vetores para armazenar valores de x0 e capacitâncias
x0_valores = linspace(x0_min, x0_max, num_pontos);
capacitancias = zeros(1, num_pontos);

% Cálculo da capacitância para cada valor de x0
for i = 1:num_pontos
    [C, ~] = Calcap(x0_valores(i));
    capacitancias(i) = C;
end

% Plotando os resultados
figure;
plot(x0_valores, capacitancias);
xlabel('Deslocamento x_0 (m)');
ylabel('Capacitância (F)');
title('Capacitância em função do deslocamento x_0');
grid on;
%% 
%Agora o cálculo da distribuição de carga com o x0= 0 e x0=5cm
Caldis(0);
Caldis(0.05);
%% 
% Potencial e campo elétrico para x0=0
epsilon0 = 8.8541e-12; % Permissividade do vácuo
ladoPlaca = 0.1; % Lado da placa (10 cm)

% Calculando a distribuição de carga
x0 = 0; % Exemplo de deslocamento de 2 cm
[rho_inferior, rho_superior] = Caldis(x0);

% Número de pontos em cada dimensão da placa
M = sqrt(length(rho_inferior));

% Redimensionando os vetores de rho para matrizes 2D
rho_matriz_inferior = reshape(rho_inferior, M, M);
rho_matriz_superior = reshape(rho_superior, M, M);

% Criando uma grade para cálculo do potencial
[X, Y] = meshgrid(linspace(0, ladoPlaca, M), linspace(0, ladoPlaca, M));
Z = 0.01; % Altura fixa acima das placas para cálculo do potencial

% Inicializando o potencial elétrico
V = zeros(size(X));

% Calculando o potencial elétrico
for i = 1:M
    for j = 1:M
        r_inferior = sqrt((X - X(i, j)).^2 + (Y - Y(i, j)).^2 + Z^2);
        V = V + rho_matriz_inferior(i, j) ./ (4 * pi * epsilon0 * r_inferior);
        r_superior = sqrt((X - (X(i, j) + x0)).^2 + (Y - Y(i, j)).^2 + Z^2);
        V = V + rho_matriz_superior(i, j) ./ (4 * pi * epsilon0 * r_superior);
    end
end
% Calculando o campo elétrico a partir do potencial
[Ex, Ey] = gradient(-V);

% Plotando o campo elétrico
figure;
quiver(X, Y, Ex, Ey);
title('Campo Elétrico Acima das Placas x_0=0');
xlabel('x (m)');
ylabel('y (m)');
axis tight;
% Plotando o potencial elétrico
figure;
surf(X, Y, V);
title('Potencial Elétrico Acima das Placas x_0=0');
xlabel('x (m)');
ylabel('y (m)');
zlabel('Potencial Elétrico (V)');
colorbar; % Adiciona uma barra de cores para representar a variação do potencial
shading interp; % Para uma visualização mais suave
%% 
% Potencial e campo elétrico x0=5cm
epsilon0 = 8.8541e-12; % Permissividade do vácuo
ladoPlaca = 0.1; % Lado da placa (10 cm)

% Calculando a distribuição de carga
x0 = 0.05; % Exemplo de deslocamento de 2 cm
[rho_inferior, rho_superior] = Caldis(x0);

% Número de pontos em cada dimensão da placa
M = sqrt(length(rho_inferior));

% Redimensionando os vetores de rho para matrizes 2D
rho_matriz_inferior = reshape(rho_inferior, M, M);
rho_matriz_superior = reshape(rho_superior, M, M);

% Criando uma grade para cálculo do potencial
[X, Y] = meshgrid(linspace(0, ladoPlaca, M), linspace(0, ladoPlaca, M));
Z = 0.01; % Altura fixa acima das placas para cálculo do potencial

% Inicializando o potencial elétrico
V = zeros(size(X));

% Calculando o potencial elétrico
for i = 1:M
    for j = 1:M
        r_inferior = sqrt((X - X(i, j)).^2 + (Y - Y(i, j)).^2 + Z^2);
        V = V + rho_matriz_inferior(i, j) ./ (4 * pi * epsilon0 * r_inferior);
        r_superior = sqrt((X - (X(i, j) + x0)).^2 + (Y - Y(i, j)).^2 + Z^2);
        V = V + rho_matriz_superior(i, j) ./ (4 * pi * epsilon0 * r_superior);
    end
end
% Calculando o campo elétrico a partir do potencial
[Ex, Ey] = gradient(-V);

% Plotando o campo elétrico
figure;
quiver(X, Y, Ex, Ey);
title('Campo Elétrico Acima das Placas x_0=5cm');
xlabel('x (m)');
ylabel('y (m)');
axis tight;
% Plotando o potencial elétrico
figure;
surf(X, Y, V);
title('Potencial Elétrico Acima das Placas x_0=5cm');
xlabel('x (m)');
ylabel('y (m)');
zlabel('Potencial Elétrico (V)');
colorbar; % Adiciona uma barra de cores para representar a variação do potencial
shading interp; % Para uma visualização mais suave
%% 
% Potencial e campo eletrico perpendicular x0= 5cm

% Parâmetros do problema
epsilon0 = 8.8541e-12; % Permissividade do vácuo
L = 200; % Dimensão da grade
x0 = 5 / 100; % Deslocamento x0 em metros (5 cm)

% Inicializando a matriz de potencial
U = zeros(L, L);

% Calculando a densidade de carga rho
[rho_placa1, rho_placa2] = Caldis(x0);

% A densidade de carga é convertida em potencial
potencial_placa1 = 220; % Potencial arbitrário para placa 1
potencial_placa2 = -220; % Potencial arbitrário para placa 2

% Posicionamento das placas na matriz de potencial
placa1_y = L/2 - 40;
placa2_y = L/2 + 40;
placa_x_start = L/2 - 80;
placa_x_end = L/2 + 80;

% Convertendo x0 para índice na matriz
x0_index = round(x0 * L / (0.1 * 2)); % Supondo que o tamanho total é 0.1*2 metros

% Aplicando potencial às placas
U(placa1_y, placa_x_start:placa_x_end) = potencial_placa1;
U(placa2_y, (placa_x_start + x0_index):(placa_x_end + x0_index)) = potencial_placa2;

% Método de relaxamento para calcular o potencial no espaço
for iter = 1:1000
    for i = 2:L-1
        for j = 2:L-1
            if U(i, j) == potencial_placa1 || U(i, j) == potencial_placa2
                continue; % Mantendo o potencial nas placas
            end
            U(i, j) = (U(i+1, j) + U(i-1, j) + U(i, j+1) + U(i, j-1)) / 4;
        end
    end
end

% Plotando o potencial elétrico
figure;
surf(U);
shading interp;
colorbar;
xlabel('Lx');
ylabel('Ly');
zlabel('Potencial Elétrico, em Volts');
title('Distribuição de Potencial Elétrico x_0=5');

% Calculando e plotando o campo elétrico
[Ex, Ey] = gradient(-U);
figure;
contour(U, 'LineWidth', 2);
hold on;
quiver(Ex, Ey, 4);
hold off;
title('Campo Elétrico x_0=5');
%% 
% Potencial e campo elétrico x0=0
epsilon0 = 8.8541e-12; % Permissividade do vácuo
L = 200; % Dimensão da grade
x0 = 0; % Deslocamento x0 em metros 

% Inicializando a matriz de potencial
U = zeros(L, L);

% Calculando a densidade de carga rho
[rho_placa1, rho_placa2] = Caldis(x0);

% A densidade de carga é convertida em potencial
potencial_placa1 = 220; % Potencial arbitrário para placa 1
potencial_placa2 = -220; % Potencial arbitrário para placa 2

% Posicionamento das placas na matriz de potencial
placa1_y = L/2 - 40;
placa2_y = L/2 + 40;
placa_x_start = L/2 - 80;
placa_x_end = L/2 + 80;

% Convertendo x0 para índice na matriz
x0_index = round(x0 * L / (0.1 * 2)); % Supondo que o tamanho total é 0.1*2 metros

% Aplicando potencial às placas
U(placa1_y, placa_x_start:placa_x_end) = potencial_placa1;
U(placa2_y, (placa_x_start + x0_index):(placa_x_end + x0_index)) = potencial_placa2;

% Método de relaxamento para calcular o potencial no espaço
for iter = 1:1000
    for i = 2:L-1
        for j = 2:L-1
            if U(i, j) == potencial_placa1 || U(i, j) == potencial_placa2
                continue; % Mantendo o potencial nas placas
            end
            U(i, j) = (U(i+1, j) + U(i-1, j) + U(i, j+1) + U(i, j-1)) / 4;
        end
    end
end

% Plotando o potencial elétrico
figure;
surf(U);
shading interp;
colorbar;
xlabel('Lx');
ylabel('Ly');
zlabel('Potencial Elétrico, em Volts');
title('Distribuição de Potencial Elétrico x_0=0');

% Calculando e plotando o campo elétrico
[Ex, Ey] = gradient(-U);
figure;
contour(U, 'LineWidth', 2);
hold on;
quiver(Ex, Ey, 4);
hold off;
title('Campo Elétrico x_0=0');
%% 

