function [Ex, Ey] = calCE(RHO, x0, x_range, y_range, dx, dy)
    % Calcula o potencial elétrico V
    V = calWi(RHO, x0, x_range, y_range, dx, dy);

    % Prepara a grade
    [X, Y] = meshgrid(min(x_range):dx:max(x_range), min(y_range):dy:max(y_range));

    % Inicializa as componentes do campo elétrico
    Ex = zeros(size(V));
    Ey = zeros(size(V));

    % Calcula o campo elétrico usando diferenças finitas
    for i = 1:size(V, 1) - 1
        for j = 1:size(V, 2) - 1
            Ex(i, j) = -(V(i, j+1) - V(i, j)) / dx;
            Ey(i, j) = -(V(i+1, j) - V(i, j)) / dy;
        end
    end

    % Plotando o campo elétrico
    figure;
    quiver(X, Y, Ex, Ey);
    title('Campo Elétrico');
    xlabel('x');
    ylabel('y');
end
