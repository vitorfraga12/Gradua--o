function [rho_placa1, rho_placa2] = Caldis(x0)
    % Constantes
    ER = 1.0;
    EO = 8.8541e-12;
    ladoPlaca = 0.1; % 10 cm
    distanciaPlacas = 0.5; % 7 cm
    N = 400;
    NT = 2 * N;
    M = sqrt(N);
    DX = ladoPlaca / M;
    DY = ladoPlaca / M;
    DL = DX;

    % Parâmetros do buraco
    ladoBuraco = 0.045; % 
    % Centralizando o buraco na placa
    distanciaVertice = 2.5/100;

    % Inicialização das matrizes e vetores
    A = zeros(NT, NT);
    B = zeros(NT, 1);
    X = zeros(NT, 1);
    Y = zeros(NT, 1);
    Z = zeros(NT, 1);

    % Preenchimento de X, Y, Z com deslocamento na placa superior
    K = 0;
    for K1 = 1:2
        for K2 = 1:M
            for K3 = 1:M
                K = K + 1;
                if K1 == 1
                    % Placa inferior
                    X(K) = DX * (K2 - 0.5);
                    Y(K) = DY * (K3 - 0.5);
                else
                    % Placa superior com deslocamento x0
                    X(K) = DX * (K2 - 0.5) + x0;
                    Y(K) = DY * (K3 - 0.5);
                end
            end
        end
    end
    Z(1:N) = 0.0;
    Z(N+1:end) = distanciaPlacas;

    % Cálculo da matriz A e vetor B
    for I = 1:NT
        for J = 1:NT
            if I == J
                A(I, J) = DL * 0.8814 / (pi * EO);
            else
                R = sqrt((X(I) - X(J))^2 + (Y(I) - Y(J))^2 + (Z(I) - Z(J))^2);
                A(I, J) = DL^2 / (4 * pi * EO * R);
            end
        end
        % Verificando se o ponto está dentro do buraco
        if X(I) > distanciaVertice && X(I) < distanciaVertice + ladoBuraco && Y(I) > distanciaVertice && Y(I) < distanciaVertice + ladoBuraco
            B(I) = 0;
        else
            B(I) = I <= N; % Ajuste na atribuição de valores para B
        end
    end

    % Calculando ρ e Capacitância
    F = inv(A);
    RHO = F * B;
    rho_placa1 = RHO(1:N); % Densidade de carga na placa inferior
    rho_placa2 = RHO(N+1:end); % Densidade de carga na placa superior
% Preparando dados para plotagem
    [X_grid, Y_grid] = meshgrid(linspace(0, ladoPlaca, M), linspace(0, ladoPlaca, M));
    rho_plot1 = reshape(rho_placa1, [M, M]);
    rho_plot2 = reshape(rho_placa2, [M, M]);

    % Plotagem da distribuição de carga
    figure;
    subplot(1, 2, 1);
    surf(X_grid, Y_grid, rho_plot1);
    title('Distribuição de Carga na Placa Inferior');
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Densidade de Carga (C/m^2)');

    subplot(1, 2, 2);
    surf(X_grid, Y_grid, rho_plot2);
    title('Distribuição de Carga na Placa Superior');
    xlabel('X (m)');
    ylabel('Y (m)');
    zlabel('Densidade de Carga (C/m^2)');

    % Retornando dados para plotagem externa, se necessário
    X_plot = X_grid;
    Y_plot = Y_grid;
    rho_plot1 = rho_plot1;
    rho_plot2 = rho_plot2;
end
