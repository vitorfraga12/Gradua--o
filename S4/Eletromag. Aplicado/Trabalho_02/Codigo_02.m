Nx=100; % Número de pontos na direção x
Ny=100; % Número de pontos na direção y
Ngrid=Nx*Ny;
dx=1; % resolução da grade
dy=1;

V=zeros(Nx,Ny); % Potencial inicial em todos os pontos

% Definindo as condições de contorno
V(1,:)=0;    
V(:,1)=0;    
V(:,end)=0;   
V(end,:)=-10;  

% Definindo os pontos da linha de potencial
x1 = 20; y1 = 87.5;
x2 = 40; y2 = 87.5;
numPoints = max(abs(x2-x1), abs(y2-y1)); % Número de pontos na linha
xLine = round(linspace(x1, x2, numPoints));
yLine = round(linspace(y1, y2, numPoints));

V=V'; % Transposição para trabalhar com a indexação correta
Vi=V(:); % Condição inicial em forma de vetor

% Formando a Matriz Laplaciana
neg_one=-1.*ones(Ngrid,1);
pos_one=1.*ones(Ngrid,1);
Lap=spdiags([pos_one pos_one 4.*neg_one pos_one pos_one], [-Nx -1 0 1 Nx], Ngrid, Ngrid);

% Aplicando as condições de contorno de Dirichlet
C=spdiags(Lap,1);
for i=1:(Nx-1)
   C(i*Nx+1)=0;
end
Lap=spdiags(C,1,Lap);
C=spdiags(Lap,-1);
for i=1:(Nx-1)
   C(i*Nx)=0;
end
Lap=spdiags(C,-1,Lap); % Matriz Laplaciana

% Definindo a linha de potencial
for i = 1:numPoints
   idx = (yLine(i)-1)*Nx + xLine(i); % Índice linear correspondente
   Lap(idx, :) = 0;  % Zera a linha na matriz Laplaciana
   Lap(idx, idx) = 1; % Define o coeficiente diagonal para 1
   Vi(idx) = 6;      % Define o valor do potencial (30V)
end

% Resolver o sistema linear
Vf=(Lap\Vi); % Potencial do capacitor
potencial=reshape(Vf,Nx,Ny)'; % conversão do vetor de potencial para matriz

%Plotando o Resultado Correto
surf(potencial);shading interp; colorbar;
xlabel('Eixo X (mm)');
ylabel('Eixo Y (mm)');
zlabel('Potencial (V)');
title('Distribuição de Potencial na Caixa');
%% 
[Ey, Ex] = gradient(-potencial, dx, dy); % O MATLAB retorna Ey antes de Ex

figure('Color','w');
quiver(Ex, Ey, 'AutoScaleFactor', 4); % Aumentar o tamanho dos vetores
axis equal tight;
title('Campo Elétrico dentro da caixa','fontsize',14);
xlabel('X','fontsize',14);
ylabel('Y','fontsize',14);
set(gca, 'FontSize', 14);
%% 

% Constante de permissividade do vácuo
epsilon0 = 8.854e-12;

% Calculando o campo elétrico normal nas arestas
% Note que você precisa extrair a componente normal do campo elétrico
% nas arestas e na linha interna. O método exato dependerá da orientação
% de cada aresta e da linha.

% Para as arestas verticais (esquerda e direita)
rho_s_left = epsilon0 * abs(Ex(:,1));
rho_s_right = epsilon0 * abs(Ex(:,end));

% Para as arestas horizontais (superior e inferior)
rho_s_top = epsilon0 * abs(Ey(1,:));
rho_s_bottom = epsilon0 * abs(Ey(end,:));

% Para a linha de potencial interna
% Você precisa calcular a componente normal do campo elétrico
% ao longo dessa linha. O cálculo exato dependerá da orientação da linha.
rho_s_line = epsilon0 * abs(Ex(20:40,70));

% Agora você tem a distribuição de carga nas arestas e na linha interna
% Você pode visualizar ou analisar esses dados conforme necessário


% Crie uma matriz de zeros para a densidade de carga
rho = zeros(Ny, Nx);

% Atribua os valores calculados de rho_s às bordas correspondentes
rho(1,:) = rho_s_top;       % Aresta superior
rho(end,:) = rho_s_bottom;  % Aresta inferior
rho(:,1) = rho_s_left;      % Aresta esquerda
rho(:,end) = rho_s_right;   % Aresta direita

% Agora atribua os valores de rho_s_line à linha interna
for i = 1:numPoints
    rho(yLine(i), xLine(i)) = rho_s_line(i);
end

% Plotando a densidade de carga superficial com surf
figure;
surf(rho);
shading interp; % Para um visual mais suave sem linhas de grade
title('Distribuição de Carga Superficial');
xlabel('X');
ylabel('Y');
zlabel('Densidade de Carga (C/m^2)');
colorbar; % Adiciona uma barra de cores para
figure;
imagesc(rho);
colormap('jet');
colorbar;
title('Distribuição de Carga a caixa');
xlabel('Coordenada X');
ylabel('Coordenada Y');
%% 
epsilon0 = 8.854e-12; % Permissividade do vácuo

% Define a região de interesse
x_range = [-50, 150]; % Limites em x
y_range = [-50, 150]; % Limites em y

% Cria uma grade para calcular o potencial
[X, Y] = meshgrid(min(x_range):dx:max(x_range), min(y_range):dy:max(y_range));
V = zeros(size(X)); % Inicializa a matriz de potencial

% Cálculo do potencial a partir da distribuição de carga
for i = 1:Nx
    for j = 1:Ny
        r = sqrt((X - i*dx).^2 + (Y - j*dy).^2);
        V = V + rho(i, j) ./ (4 * pi * epsilon0 * r);
    end
end

% Tratar singularidades
V(isinf(V)) = NaN;

% Calculando o campo elétrico
[Ex, Ey] = gradient(-V, dx, dy);

% Plotando o potencial
figure;
surf(X, Y, V);
shading interp;
colorbar;
title('Distribuição de Potencial');
xlabel('x');
ylabel('y');
zlabel('Potencial (V)');

% Plotando o campo elétrico
figure;
quiver(X, Y, Ex, Ey, 'k-', 'LineWidth', 1.5);
title('Campo Elétrico');
xlabel('x');
ylabel('y');