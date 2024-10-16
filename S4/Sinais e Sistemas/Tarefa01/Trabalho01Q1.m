%% Questão 01
hz = 60;  % f = 60Hz
periodo = 1 / hz;  % T = 0.0166666667s
pvalores = [3, 5, 10, 50];  % Numéro de amostras para cada período

% Plotamos uma senoide para cada valor de P, sendo eles 3, 5, 10 e 50
for p = pvalores
    n = 0:(p-1);  % Vetor de índices de tempo discreto
    t = n * periodo / p;  % Vetor de tempo discreto
    
    sen = sin(2*pi*hz*t);  % Senoide no domínio do tempo discreto
    
    % Plote a senoide
    figure;
    stem(n, sen, 'k', 'filled');
    title([num2str(p) ' amostras por período']);
    xlabel('Amostras');
    ylabel('Amplitude');
    grid on;
end


