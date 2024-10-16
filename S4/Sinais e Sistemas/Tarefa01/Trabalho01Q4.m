%% Questão 04
% Parâmetros
n = -15:15;  % Amostras de tempo discreto
t = 10; % Período

% Criamos uma função impulso unitário e logo em seguida usamos ela para
% criar uma função degrau unitário
impulso = zeros(size(n));
impulso(n == 0) = 1;
degrau_impulso = cumsum(impulso);

% Criamos uma função degrau unitário e logo em seguida usamos ela para
% criar uma função impulso unitário
impulso_degrau = diff(degrau_impulso);
impulso_degrau = [0, impulso_degrau];

% Não consegui gerar uma onda quadrada usando a função de Degrau Unitário,
% então fiz ela de outra maneira
onda_quadrada = 2 * (mod(n, t) < t/2) - 1;

% Plote os gráficos
figure;

% Impulso unitário
subplot(3, 1, 1);
stem(n, impulso_degrau,'k', 'filled');
title('Função Impulso Unitário');
xlabel('Tempo discreto (n)');
ylabel('Amplitude');
grid on;

% Degrau unitário
subplot(3, 1, 2);
stem(n, degrau_impulso,'k', 'filled');
title('Função Degrau Unitário');
xlabel('Tempo discreto (n)');
ylabel('Amplitude');
grid on;

% Onda quadrada
subplot(3, 1, 3);
stem(n, onda_quadrada,'k', 'filled');
title('Onda Quadrada a partir do Degrau');
xlabel('Tempo discreto (n)');
ylabel('Amplitude');
grid on;
