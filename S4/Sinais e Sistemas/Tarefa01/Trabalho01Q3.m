%% Questão 3
% Parametros
n = 0:200;
w0 = [0.1*pi, 1*pi, 2*pi, 4*pi]; % Utilizei a multiplicação por pi pois o livro informa que 𝜔0 precisa ter unidade de radiano.

% Matriz para armazenar os diferentes valores que utilizamos para 𝜔0;
x = zeros(length(w0), length(n));

% Calculamos os sinais para os diferentes valores de 𝜔0
for i = 1:length(w0)
    x(i, :) = exp(1j * w0(i) * n);
end

% Plotamos os sinais para os diferentes valores de 𝜔0.
figure;
for i = 1:length(w0)
    subplot(length(w0), 1, i);
    stem(n, x(i, :));
    title(['𝜔0 = ', num2str(w0(i))]);
    xlabel('n');
    ylabel('Amplitude');
end
