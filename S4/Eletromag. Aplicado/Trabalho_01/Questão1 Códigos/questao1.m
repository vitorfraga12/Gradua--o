ER = 1.0;
EO = 8.8541e-12;
N = 100; % Número de segmentos em cada fio
VO = 2.0; % Diferença de potencial

% Define o intervalo para x0
x0_values = 0:0.1:5; % Exemplo: de 0 a 5 com passo de 0.1

% Inicializa as capacitâncias
C_total = zeros(size(x0_values));

% Loop para diferentes valores de x0
for idx = 1:length(x0_values)
    x0 = x0_values(idx);

    % Capacitor 1 (segmento horizontal inicial)
    d1 = 7.0; % Distância entre os fios
    l1 = 2.0+x0; % Comprimento do fio (fixo para a parte horizontal inicial)
    C1 = calCap(d1, l1, x0, N, EO);

    % Capacitor 2 (segmento horizontal final)
    d2 = 7.75; % Distância entre os fios
    l2 = 0.75; % Comprimento do fio (fixo para a parte horizontal final)
    C2 = calCap(d2, l2, x0, N, EO);

    % Capacitor 2 (segmento horizontal final)
    d3 = 7.0; % Distância entre os fios
    l3 = 7.25-x0; % Comprimento do fio (fixo para a parte horizontal final)
    C3 = calCap(d3, l3, x0, N, EO);

    % Soma as capacitâncias
    C_total(idx) = C1 + C2 + C3;
end

% Exibe a capacitância total em função de x0
plot(x0_values, C_total);
xlabel('x_0 (Deslocamento)');
ylabel('Capacitância Total');
title('Capacitância Total em Função do Deslocamento x_0');

% agora vemos com a distribuição de carga nos fios
% x0= 0.5
x0 = 0.5; % Valor de exemplo para x0
d = 7; % Distância entre os fios
l = 10.0; % Comprimento do fio
N = 100; % Maior número de segmentos para melhor resolução
EO = 8.8541e-12;

RHO = calR(d, l, x0, N, EO);

% Posições dos segmentos
DL = l / N;
X = (DL/2):DL:l;

% Plot
figure;
plot(X, RHO(1:N), 'b-', X, RHO(N+1:end), 'r--');
legend('Fio Inferior', 'Fio Superior');
xlabel('Posição (m)');
ylabel('Densidade de Carga (C/m)');
title('Distribuição de Carga no Capacitor x_0=0.45');

% x0= 5
x0 = 5; % Valor de exemplo para x0
d = 7; % Distância entre os fios
l = 10.0; % Comprimento do fio
N = 100; % Maior número de segmentos para melhor resolução
EO = 8.8541e-12;

RHO = calR(d, l, x0, N, EO);

% Posições dos segmentos
DL = l / N;
X = (DL/2):DL:l;

% Plot
figure;
plot(X, RHO(1:N), 'b-', X, RHO(N+1:end), 'r--');
legend('Fio Inferior', 'Fio Superior');
xlabel('Posição (m)');
ylabel('Densidade de Carga (C/m)');
title('Distribuição de Carga no Capacitor x_0=5');

%Agora para p na geometria 
%x0= 0.5
x0=0.5
gerageo(RHO, x0);

%x0= 5
x0=5;
gerageo(RHO, x0);
%% 
%por fim potencial e campo elétrico
%x0=0.5
x_range = [-5, 15];
y_range = [-10, 10];
dx = 0.1;
dy = 0.1;
x0= 0.5;
calWi(RHO, x0, x_range, y_range, dx, dy);
calCE(RHO, x0, x_range, y_range, dx, dy);


%% 
%x0=5
x_range = [-5, 15];
y_range = [-10, 10];
dx = 0.1;
dy = 0.1;
x0= 5;
calWi(RHO, x0, x_range, y_range, dx, dy);
calCE(RHO, x0, x_range, y_range, dx, dy);

%% 


