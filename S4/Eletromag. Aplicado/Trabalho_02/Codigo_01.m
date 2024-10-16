L = 10;
discretizacao = 100;
V1 = 10; % Para a primeira linha
V0 = 6; % De (2, 7) até (4, 7)
U = zeros(80, discretizacao);  % Matriz 80x100
U(80, 1:100) = V1;  % Potencial de 10 V na borda superior
U(70, 20:40) = V0;

% Loop para calcular o potencial
for k = 1:10000
    for i = 2:79
        for j = 2:99
            if U(i, j) == V1 || U(i, j) == V0 % Mantém o potencial nas bordas e na linha
                continue;
            end
            U(i, j) = (U(i - 1, j) + U(i + 1, j) + U(i, j - 1) + U(i, j + 1)) / 4;
        end
    end
end
surf(U);shading interp; colorbar;
xlabel('Eixo X (mm)');
ylabel('Eixo Y (mm)');
zlabel('Potencial (V)');
title('Distribuição de Potencial na Caixa');
%% 
[Ex,Ey]=gradient(-U);
epi = 8.85e-12;
% Calcula a divergência do campo elétrico
[divEx, ~] = gradient(Ex);
[~, divEy] = gradient(Ey);
divE = abs(divEx) + abs(divEy);
% Calcula a densidade de carga
rho = epi .* divE;
% Visualização da distribuição de carga
figure;
surf(rho); shading interp; colorbar;
xlabel('Eixo X (mm)');
ylabel('Eixo Y (mm)');
zlabel('\rho');
title('Distribuição de Carga na Caixa');
% Visualização alternativa da distribuição de carga
figure;
imagesc(rho); colormap('jet'); colorbar;
title('Distribuição de Carga na Caixa');
xlabel('Eixo X (mm)');
ylabel('Eixo Y (mm)');
figure, contour(U,'LineWidth',2);
hold on, quiver(Ex,Ey,4), hold off
axis tight;
xlabel('Eixo X (mm)');
ylabel('Eixo Y (mm)');
title('Campo Elétrico na Estrutura')
%% 
dx = 0.5; 
dy = 0.5;
% Define a região de interesse fora da caixa
x_f = [-20, 80]; % Limites em x
y_f = [-20, 80]; % Limites em y

% Cria uma grade para calcular o potencial fora da caixa
[X, Y] = meshgrid(min(x_f):dx:max(x_f), min(y_f):dy:max(y_f));
V_exterior = zeros(size(X)); % Inicializa a matriz de potencial exterior

% Cálculo do potencial exterior a partir da distribuição de carga
for i = 1:size(U, 1)
    for j = 1:size(U, 2)
        r = sqrt((X - i*dx).^2 + (Y - j*dy).^2);
        V_exterior = V_exterior + rho(i, j) ./ (4 * pi * epi * r);
    end
end

% Tratar singularidades
V_exterior(isinf(V_exterior)) = NaN;

% Calculando o campo elétrico exterior
[Ex_exterior, Ey_exterior] = gradient(-V_exterior, dx, dy);

% Plotando o potencial exterior
figure;
surf(X, Y, V_exterior);
shading interp;
colorbar;
title('Distribuição de Potencial Exterior da Caixa');
xlabel('Eixo X (mm)');
ylabel('Eixo Y (mm)');
zlabel('Potencial (V)');

% Plotando o campo elétrico exterior
figure;
quiver(X, Y, Ex_exterior, Ey_exterior, 'k-', 'LineWidth', 1.5);
title('Campo Elétrico Exterior');
xlabel('Eixo X (mm)');
ylabel('Eixo Y (mm)');