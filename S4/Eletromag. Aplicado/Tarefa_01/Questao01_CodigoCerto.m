% Valores de a para calcular a integral
valores_a = [2, 5];

% Números de segmentos para discretização
segmentos = 2:2:30000;

% Inicialize vetores para armazenar erros para cada valor de a
erros_a2 = zeros(size(segmentos));
erros_a5 = zeros(size(segmentos));

% Calcule o valor exato da integral
valor_exato = pi/2;

% Loop através dos valores de a
for i = 1:length(valores_a)
    a = valores_a(i);
    
    % Loop através dos números de segmentos
    for j = 1:length(segmentos)
        n = segmentos(j);
        h = a / n;
        
        % Regra de Euler
        integral_euler = 0;
        for k = 0:n-1
            x = k * h;
            integral_euler = integral_euler + h * (1/sqrt(a^2 - x^2));
        end
        
        % Regra de 1/3 Simpson
        integral_simpson = 0;
        for k = 0:n-1
            x = k * h;
            if k == 0 || k == n-1
                integral_simpson = integral_simpson + (h/3) * (1/sqrt(a^2 - x^2));
            elseif mod(k, 2) == 1
                integral_simpson = integral_simpson + (h/3) * (4 * (1/sqrt(a^2 - x^2)));
            else
                integral_simpson = integral_simpson + (h/3) * (2 * (1/sqrt(a^2 - x^2)));
            end
        end
        
        % Calcule o erro para cada método
        erro_euler = abs(integral_euler - valor_exato);
        erro_simpson = abs(integral_simpson - valor_exato);
        
        % Armazene os erros nos vetores
        if a == 2
            erros_a2(j) = erro_euler;
        elseif a == 5
            erros_a5(j) = erro_euler;
        end
    end
end

% Crie gráficos dos erros em função do número de segmentos
figure;
plot(segmentos, erros_a2, 'r-', 'LineWidth', 2, 'DisplayName', 'a = 2');
hold on;
plot(segmentos, erros_a5, 'b-', 'LineWidth', 2, 'DisplayName', 'a = 5');
xlabel('Número de Segmentos');
ylabel('Erro Absoluto');
legend('Location', 'northwest');
title('Erro da Integração Numérica');
grid on;


