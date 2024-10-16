c = 3e8; 
f = 575e6; 
lambda = c / f; 
E0 = 1; 
L = 1;
theta = 60 * pi / 180;
z = linspace(0, L, 1000);

% Calculando a força eletromotriz (fem) para cada lado inclinado
fem = zeros(size(z)); % Inicialização do vetor de fem

for i = 1:length(z)
    fem(i) = 2 * E0 * z(i) * cos(theta); % Calculando a fem para cada ponto ao longo do lado inclinado
end

% Plotagem do gráfico
figure;
plot(z, fem);
xlabel('Distância ao longo dos lados inclinados (m)');
ylabel('Força Eletromotriz (V)');
title('Força Eletromotriz Induzida');

% Encontrando a faixa onde a fem cai para 1/sqrt(2) de seu valor máximo
max_fem = max(fem);
threshold = max_fem / sqrt(2);
indices = find(fem >= threshold);
min_index = indices(1);
max_index = indices(end);

hold on;
plot(z(min_index), fem(min_index), 'ro');
plot(z(max_index), fem(max_index), 'ro');
text(z(min_index), fem(min_index), sprintf('(%.2f m, %.2f V)', z(min_index), fem(min_index)), 'VerticalAlignment', 'bottom');
text(z(max_index), fem(max_index), sprintf('(%.2f m, %.2f V)', z(max_index), fem(max_index)), 'VerticalAlignment', 'top');
hold off;
grid on;