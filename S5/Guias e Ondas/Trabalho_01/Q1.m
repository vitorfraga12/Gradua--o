% Primeira parte
% Definindo os parametros a ser usado
discretization = 1000;
c =  physconst('lightspeed'); %Váriavel da Velocidade da Luz

f = linspace(1e6, 1e9, discretization); % Intervalo de frequência
e0 = 1; % Amplitude
b = 0.17647058823; % Comprimento da antena a partir da relação b=c/2f
z1 = 1;
t_seconds = 200;

% Cálculos necessários para encontrar a f.e.m induzida 
lamb = c ./ f; % Comprimento de onda 
k = 2*pi ./ lamb; % Número de onda
frequencia_ang = 2*pi*f; %Frequência Angular
phi = (frequencia_ang .* t_seconds) - (k - z1); 

% Cálculo da f.e.m induzida
fem = (sqrt(3) * b) * e0 * (-sin(phi - (k .* b) / 2));

% Cálculo do valor máximo que a f.e.m atinge
max_fem = max(abs(fem));

% Cálculo do valor de limiar para 1/sqrt(2) do valor máximo
threshold = max_fem / sqrt(2);

% Gráfico
plot(f, fem);
hold on;
yline(threshold, 'r--', 'Faixa');
fill([f, fliplr(f)], [fem, zeros(size(fem))], 'b', 'FaceAlpha', 0.5);
hold off;
xlabel('Frequência [Hz]');
ylabel('F.E.M [V]');
title('F.E.M vs Frequência');
legend('F.E.M', 'Faixa', 'Location', 'best');
grid on
%% 
% Atualização a fase
phi = frequencia_ang .* t_seconds - (k .* z1) + pi/6;

% Recalculando todos os valores necessários
fem = (sqrt(3) * b) * e0 * (-sin(phi - (k .* b) / 2));
max_fem = max(abs(fem));
threshold = max_fem / sqrt(2);

% Achando a banda de frequência onde a f.e.m fica inferior ao limiar
below_threshold = f(fem < threshold);

% Gráfico
plot(f, fem);
hold on;
yline(threshold, 'r--', 'Faixa');
fill([f, fliplr(f)], [fem, zeros(size(fem))], 'g', 'FaceAlpha', 0.5);
hold off;
xlabel('Frequência [Hz]');
ylabel('F.E.M [V]');
title('F.E.M vs Frequência');
legend('F.E.M', 'Faixa', 'Location', 'best');
grid on

