function RHO = calR(d, l, x0, N, EO)
    NT = 2 * N;
    DL = l / N;
    A = zeros(NT, NT);
    B = ones(NT, 1);
    B(N+1:end) = -1.0;

    % Define as posições dos segmentos nos fios
    for K = 1:N
        X(K) = (K - 0.5) * DL;
        Y(K) = 0;
        Z(K) = 0;
        
        % Calcula a posição X para o fio superior considerando o "dente"
        if X(K) < 2
            X(K+N) = X(K) + x0; % Parte horizontal inicial
        elseif X(K) < 2.75
            X(K+N) = 2 + x0; % Parte vertical
        else
            X(K+N) = X(K) - 0.75 + x0; % Parte horizontal final
        end

        % Calcula a posição Y para o fio superior considerando o "dente"
        if X(K) < 2
            Y(K+N) = d; % Parte horizontal inicial
        elseif X(K) < 2.75
            Y(K+N) = d + (X(K) - 2); % Parte vertical
        else
            Y(K+N) = d + 0.75; % Parte horizontal final
        end

        Z(K+N) = 0;
    end

    % Constrói a matriz de coeficientes
    for I = 1:NT
        for J = 1:NT
            if I == J
                A(I, J) = DL * 0.8814 / (pi * EO);
            else
                R = sqrt((X(I) - X(J))^2 + (Y(I) - Y(J))^2 + (Z(I) - Z(J))^2);
                A(I, J) = DL^2 / (4 * pi * EO * R);
            end
        end
    end

    % Calcula a distribuição de carga e a capacitância
    F = inv(A);
    RHO = F * B;
end