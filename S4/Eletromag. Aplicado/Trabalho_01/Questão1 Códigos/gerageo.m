function gerageo(RHO, x0)
    % Define os parâmetros
    L = 10;
    d = 7;
    L2 = 2;
    N = length(RHO) / 2; % Número de segmentos em cada fio
    DL = L / N;

    % Cria a figura
    figure;
    hold on;
    axis equal;
    grid on;
    xlabel('x');
    ylabel('y');

    % Desenha o fio inferior com carga
    for i = 1:N
        x = (i - 0.5) * DL;
        scatter(x, 0, 10, RHO(i), 'filled'); % Carga negativa no fio inferior
    end

    % Desenha o fio superior com "dente" e carga
    for i = 1:N
        x = (i - 0.5) * DL;

        % Calcula a posição do fio superior considerando o "dente"
         if x < 2
            xUpper = x + x0; % Parte horizontal inicial
            yUpper = d; % Parte horizontal inicial
        elseif x < 2.75
            xUpper = 2 + x0; % Parte vertical
            yUpper = d + (x - 2); % Parte vertical
         elseif x<3.5
            xUpper = x - 0.75 + x0; % Parte horizontal final
            yUpper = d + 0.75; % Parte horizontal final
         elseif x < 4.25
            xUpper = 2.75 + x0; % Parte vertical
            yUpper = 1.6*d -x; % Parte vertical
         else
             xUpper = x-1.5 +x0; % Parte horizontal inicial
             yUpper = d; % Parte horizontal inicial
         end

        scatter(xUpper, yUpper, 10, RHO(i+N), 'filled'); % Carga positiva no fio superior
    end

    % Limites do gráfico e ajustes finais
    colormap jet; % Define uma escala de cores
    colorbar; % Adiciona uma barra de cores para referência
    xlim([-1, L+1]);
    ylim([-1, d+2]);

    hold off;
end
