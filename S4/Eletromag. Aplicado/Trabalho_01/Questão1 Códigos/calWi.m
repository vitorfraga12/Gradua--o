function V = calWi(RHO, x0, x_range, y_range, dx, dy)
    % Constantes
    epsilon0 = 8.854e-12; % Permissividade do vácuo
    L = 10; % Comprimento total dos fios
    d = 7; % Distância vertical entre os fios
    N = length(RHO) / 2; % Número de cargas em cada fio
    DL = L / N; % Comprimento de cada segmento

    % Cria uma grade para calcular o potencial
    [X, Y] = meshgrid(min(x_range):dx:max(x_range), min(y_range):dy:max(y_range));
    V = zeros(size(X)); % Inicializa a matriz de potencial

    % Loop sobre as cargas no fio inferior
    for i = 1:N
        x_charge = (i - 0.5) * DL;
        y_charge = 0;
        r = sqrt((X - x_charge).^2 + (Y - y_charge).^2);
        V = V + RHO(i) ./ (4 * pi * epsilon0 * r);
    end

    % Loop sobre as cargas no fio superior (com a geometria específica)
    for i = 1:N
        x = (i - 0.5) * DL;
        
        % Define a posição y do fio superior, com base na geometria dada
        if x < 2
            xUpper = x + x0; yUpper = d;
        elseif x < 2.75
            xUpper = 2 + x0; yUpper = d + (x - 2);
        elseif x < 3.5
            xUpper = x - 0.75 + x0; yUpper = d + 0.75;
        elseif x < 4.25
            xUpper = 2.75 + x0; yUpper = 1.6 * d - x;
        else
            xUpper = x - 1.5 + x0; yUpper = d;
        end

        r = sqrt((X - xUpper).^2 + (Y - yUpper).^2);
        V = V + RHO(i + N) ./ (4 * pi * epsilon0 * r);
    end

    % Lidar com singularidades
    V(isinf(V)) = NaN;

    % Plotando o potencial
    figure;
    surf(X, Y, V);
    shading interp;
    colorbar;
    title('Distribuição de Potencial');
    xlabel('x');
    ylabel('y');
    zlabel('Potencial (V)');
end