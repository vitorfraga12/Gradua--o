% Parâmetros
x = 0:0.1:10; % Valores de x no intervalo [0, 10]
m_values = [0, 1, 2]; % Valores de m

% Inicialização dos resultados
Jm_approx = zeros(length(x), length(m_values));

% Cálculo das funções de Bessel usando o método de 1/3 Simpson
for i = 1:length(m_values)
    m = m_values(i);
    for j = 1:length(x)
        integral_result = 0;
        N = 1000; % Número de pontos para a discretização
        a = 0;
        b = pi;
        h = (b - a) / N;
        
        for k = 0:N
            theta = a + k * h;
            if k == 0 || k == N
                coefficient = 1;
            elseif mod(k, 2) == 1
                coefficient = 4;
            else
                coefficient = 2;
            end
            integral_result = integral_result + coefficient * cos(x(j) * sin(theta) - m * theta);
        end
        
        integral_result = integral_result * h / 3;
        Jm_approx(j, i) = integral_result / pi;
    end
end

% Cálculo das funções de Bessel usando a função do MATLAB
Jm_exact = zeros(length(x), length(m_values));
for i = 1:length(m_values)
    m = m_values(i);
    for j = 1:length(x)
        Jm_exact(j, i) = besselj(m, x(j));
    end
end

% Plot dos resultados
figure;

for i = 1:length(m_values)
    subplot(length(m_values), 1, i);
    plot(x, Jm_approx(:, i), 'r-', x, Jm_exact(:, i), 'g--');
    title(['J' num2str(m_values(i)) '(x) Aproximado vs. Exato']);
    legend('Aproximado', 'Exato');
    grid on;
end

% Cálculo e plot do erro
error = abs(Jm_exact - Jm_approx);
figure;
plot(x, error, 'LineWidth');
title('Erro de Aproximação');
legend('Erro J0', 'Erro J1', 'Erro J2');
grid on;