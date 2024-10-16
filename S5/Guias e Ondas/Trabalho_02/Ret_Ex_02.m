% Parâmetros do sinal modulado AM
A_c = 1; % Amplitude da portadora
k = 0.5; % Índice de modulação
fs = 1e9; % Frequência de amostragem (Hz)
t = 0:1/fs:1e-5; % Tempo (s)
f_carr = 10e9; % Frequência da portadora (Hz)
f_mod = 100e6; % Frequência de modulação (Hz)

% Sinal de mensagem m(t)
m_t = sin(2 * pi * f_mod * t);

% Sinal modulado em AM com a equação s(t) = Ac[1 + k m(t)] cos(2πf_c t)
s_t = A_c * (1 + k * m_t) .* cos(2 * pi * f_carr * t);

% Plote o sinal modulado
figure;
plot(t(1:500), s_t(1:500));
title('Sinal Modulado em Amplitude (AM)');
ylim([0,2])
xlabel('Tempo (s)');
ylabel('Amplitude');

%% Parâmetros do guia de onda
l = 0.025; % Largura do guia de onda (m)
h = 0.01; % Altura do guia de onda (m)
c = 3e8;    % Velocidade da luz (m/s)
f_op = 10e9; % Frequência de operação (Hz)

% Definir o modo TE10 e calcular a frequência de corte
m = 1;
n = 0;
fc = (c / 2) * sqrt((m / l)^2 + (n / h)^2);

% Verificar se o modo é propagante
if f_op > fc
    fprintf('O modo TE10 é propagante.\n');
else
    fprintf('O modo TE10 não é propagante.\n');
end

% Simulação da atenuação do sinal ao longo do guia
attenuation_factor = exp(-t .* 5e6); % Atenuação mais lenta
propagated_signal = s_t .* attenuation_factor;

% Plotar o sinal propagado
figure;
plot(t(1:500), propagated_signal(1:500));
title('Sinal Propagado no Guia de Onda');
ylim([0,2])
xlabel('Tempo (s)');
ylabel('Amplitude');