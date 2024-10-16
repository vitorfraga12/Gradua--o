% Parâmetros
n = 0:200;  % Amostras de tempo usada durante todo o código
c = 2 + 1i;  % Escolha de fixa de C para todo o código

% Parâmetros diferente para a variavel "a" para cada questão
a_a = -0.01 + 0i;  % Escolha para o item a). A parte imaginária precisa ser zerada.
a_b = 0 + 0.1i; % Escolha para o item b). A parte real precisa ser zerada.
a_c = 0.01 + 0.1i; % Escolha para o item c). As duas partes não podem estar zeradas.

% Sinais que usaremos para cada item.
x_a = c * exp(a_a * n);
x_b = c * exp(a_b * n);
x_c = c * exp(a_c * n);
%% Item a)

stem(n, x_a);
xlabel('Numéro de amostras');
ylabel('Amplitude do Sinal');
title('Parte Real do Sinal Exponencial');

%% Item b)

stem(n, x_b);
xlabel('Numéro de amostras');
ylabel('Amplitude do Sinal');
title('Sinal Oscilatório');

%% Item c)

stem(n, x_c);
xlabel('Numéro de amostras');
ylabel('Amplitude do Sinal');
title('Sinal Amortecido');
