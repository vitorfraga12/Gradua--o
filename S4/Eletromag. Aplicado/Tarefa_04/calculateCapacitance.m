function C = calculateCapacitance(d, l, x0, N, EO)
    NT = 2 * N;
    DL = l / N;
    A = zeros(NT, NT);
    B = ones(NT, 1);
    B(N+1:end) = -1.0;

    % Define as posições dos segmentos nos fios
    for K = 1:N
        X(K) = (K - 0.5) * DL;
        X(K+N) = (K - 0.5) * DL + x0; % Deslocamento aplicado ao fio superior
        Y(K) = 0;
        Y(K+N) = d;
        Z(K) = 0;
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
    Q = sum(RHO(1:N)) * DL;
    C = abs(Q) / 2;
end

