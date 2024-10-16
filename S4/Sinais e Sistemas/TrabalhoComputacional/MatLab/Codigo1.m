%% Começamos o código lendo e plotando o sinal de áudio dado pelo professor
[x, fs] = audioread('ImperialPlusCantina.wav');
plot(x);
xlabel('Amostras');
title('Sinal de Áudio');
%% Em seguida deixamos o áudio no dominio da frequência e plotamos sua DEP
x_n = x;    
X_ejw = fft(x_n);
X_ejw_s = fftshift(X_ejw);

% Prepara o eixo de frequência centralizado
tamanho = length(x_n);              % Número de pontos no sinal
f_geral = (-tamanho/2:tamanho/2-1)*(fs/tamanho); % Eixo de frequência centralizado usaremos f_geral pois é o mesmo tamanho para todos os sinais

% Calcula a magnitude e ajusta para a parte centralizada
magnitude_x = abs(X_ejw_s)/max(real(X_ejw_s));

% Plota o espectro centralizado
plot(f_geral, magnitude_x);
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Densidade Espectral de Potência do Áudio');
%% Para fazer o filtro passa-baixa é necessario que h[n] seja uma função sinc, então declamos a função h[n]
% Temos tanto fs quanto o número de amostras
% Calculamos a duração do sinal
duracao = (tamanho-1) / fs;

% Criamos n em torno de 0, para obtermos a função sinc
n = (-duracao/2 : 1/fs : duracao/2);

% Frequência de corte de 4000Hz
fc = 4000; 

% Calculando a função sinc
h_n = sinc( 2 * fc * n);
% Tornando h_n um vetor coluna para conseguirmos multiplicar e
% convolucionar com X_ejw e x_n respectivamente
h_n = h_n(:);

% Plotando a função sinc
plot(n, h_n);
xlabel('Amostras');
title('Sinal h[n]');
%% Agora iremos deixar a função h[n] no domínio da frequência
H_ejw = fft(h_n);

% Centralizando o espectro
H_ejw_s = fftshift(H_ejw);

% Calculando a magnitude
magnitude_h = abs(H_ejw_s)/max(real(H_ejw_s));

% Plotando o espectro da função sinc
plot(f_geral, magnitude_h);
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Densidade Espectral de Potência H(e^{j\omega})');
%% Multiplicamos X(e^jw) por H(e^jw) para conseguir obter o sinal filtrado e logo em seguida plotamos a DEP do sinal
Y_ejw = (X_ejw) .* (H_ejw);
Y_ejw_s = fftshift(Y_ejw);
magnitude_y = abs(Y_ejw_s)/max(real((Y_ejw_s)));
plot(f_geral, magnitude_y);
xlabel('Frequência (Hz)');
ylabel('Magnitude');
title('Densidade Espectral de Potência Y(e^{j\omega})');
%% Transformamos do dominio da frequência para o domínio do tempo, para que assim tenhamos o áudio
% Calculando a IFFT para obter o sinal filtrado no domínio do tempo
y_n1 = ifft(Y_ejw);
%Quando eu faço isso, obtenho um sinal onde a primeira metade do áudio
%começa no final e a segunda metade começa no inicio, então altero o valor
%de y_n1 para que eu possa ficar com o áudio no tempo correto

tamy = length(y_n1); % Obter o tamanho do vetor
y_n1 = [y_n1((tamy/2)+1:end); y_n1(1:tamy/2)];

% Plotando o sinal filtrado
plot(real(y_n1)); % Usar 'real' para garantir que qualquer parte imaginária residual seja descartada
xlabel('Amostras');
title('Sinal com Filtragem na Frequência');
audiowrite('Filtro_Multiplicacao.wav', y_n1, fs)
%% Para esse caso, temos que a multiplicação no dominio da frequência é igual a convolução no dominio do tempo, então fazemos a convolução de x[n] com h[n]
y_n2 = conv(x_n, h_n);
%Quando fazemos a convolução completa temos que os 30s finais e iniciais
%não tem nada por conta da sobreposição, com isso vou limitar y_n2 apenas
%para a parte com contéudo
elements_30_sec = 30 * fs;
y_n2 = y_n2(elements_30_sec + 1 : end - elements_30_sec);
plot(y_n2);
xlabel('Amostras');
title('Sinal com Filtragem no Tempo');

% Ouvindo o sinal resultante
audiowrite('Filtro_Convolucao.wav',y_n2, fs);
%% Vamos utilizar da função que cria o filtro passa baixa para compararmos com o nosso filtro
frequencia_de_corte = 4000; 

% Criar um filtro passa-baixa   
filter = designfilt('lowpassfir', 'FilterOrder', 20, ...
               'CutoffFrequency', frequencia_de_corte, 'SampleRate', fs);

% Aplicar o filtro ao sinal
sinal_filtrado = filtfilt(filter, x);

% Salvar o sinal filtrado
audiowrite('Filtro_FuncaoPassaBaixa.wav', sinal_filtrado, fs);