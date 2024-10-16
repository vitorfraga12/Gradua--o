% Parâmetros do guia de onda retangular
a = 5e-2; % Largura do guia (m)
b = 2e-2; % Altura do guia (m)
c = 3e8; % Velocidade da luz no vácuo (m/s)

% Número máximo de modos para cálculo
max_m = 5;
max_n = 5;

% Inicialização da matriz para armazenar as frequências de corte
frequencias_corte = zeros(max_m, max_n);

% Cálculo das frequências de corte para os modos TE_{mn}
for m = 0:max_m
    for n = 0:max_n
        if m == 0 && n == 0
            frequencias_corte(m+1, n+1) = NaN; % Modo TEM não é suportado
        else
            frequencias_corte(m+1, n+1) = (c/2) * sqrt((m/a)^2 + (n/b)^2);
        end
    end
end

% Plot das frequências de corte
figure;
imagesc(frequencias_corte);
colorbar;
title('Dispersão das Frequências de Corte para Modos TE_{mn}');
xlabel('Número de Modos na Direção n');
ylabel('Número de Modos na Direção m');
set(gca, 'XTick', 1:max_n, 'XTickLabel', 0:max_n);
set(gca, 'YTick', 1:max_m, 'YTickLabel', 0:max_m);

% Mostrar os valores das frequências de corte na tabela
for m = 0:max_m
    for n = 0:max_n
        text(n+1, m+1, sprintf('%.2f GHz', frequencias_corte(m+1, n+1)*1e-9), ...
            'HorizontalAlignment', 'center', 'Color', 'w');
    end
end
