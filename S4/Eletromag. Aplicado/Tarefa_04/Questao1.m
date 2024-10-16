L = 10;
discretizacao = 100;
U = zeros(discretizacao, discretizacao);
line_size = 1;
U(100, 1:100) = 10;
count = 1;
tolerance = 0.3;

% loops
for k = 1:10000
    for i = 2:discretizacao - 1
        for j = 2:discretizacao - 1
            % Let the voltage's conditions intact
            if U(i, j) == 1 || U(i, j) == -1
                continue;
            end
            % Laplace equation: Delta U = 0
            U(i, j) = (U(i - 1, j) + U(i + 1, j) + U(i, j - 1) + U(i, j + 1)) / 4;
        end
    end
end
% plot the results
surf(U);shading interp; colorbar;
xlabel('Eixo X');
ylabel('Eixo Y');
zlabel('Potencial (V)');
title ('Potencial Elétrico no Capacitor de Placas Paralelas');
view(-54,6);
%  Electrostatic field
[Ex,Ey]=gradient(U);
epi=8.85e-12;
rho=epi*(Ex^2+Ey^2);
figure, contour(U,'LineWidth',2);
hold on, quiver(Ex,Ey,4), hold off
title('Campo elétrico do Capacitor em caso estático');
%Thrid figure for better viusalization of the electric field
figure,contour(U,'LineWidth',2); hold on, quiver(Ex,Ey,4), hold off;
zoom(3);
axis([10 40 130 160]);
xlabel('Eixo X'), ylabel('Eixo Y'), title(' Electrcial field :Zoom x4 Left Corner of the 1st plate ')
% Plote a distribuição de carga
figure;
surf(rho);shading interp; colorbar;
xlabel('Eixo X');
ylabel('Eixo Y');
zlabel('Potencial (V)');
title ('Potencial Elétrico do Capacitor');
figure;
imagesc(rho);
colormap('jet');
colorbar;
title('Distribuição de Carga no Fio');
xlabel('Eixo X');
ylabel('Eixo Y');
%% 
syms m n
L = 1;                  % Length of wire
lim1 = 0;               % Lower limit
lim2 = L;               % Upper limit
a = 0.001;              % radius of the wire
V = 10;                  % Potential V = 1 volt
ep_0 = 8.854*1e-12;     % permittivity
ep_r = 1;               % relative Permittivity
N = 10;                 % Number of sections
M = N;                  % Number of testing function
hx = (lim2-lim1)/(N);   % Step width of each section
x = lim1:hx:lim2;
dx = hx;
% Definition of testing function
Zmn = zeros(M,N);
Zmn1 = zeros(M,N);
%Basis function
% 1 for (n-1)delx < x < n*delx
% fn = 1;
ii = 1:N;
xm = (ii-0.5)*dx;        % Center for delta function
for n = 1:N
    for m = 1:N
        
        xb = n*dx;              % Section higher limit
        xa = (n-1)*dx;          % Section lower limit
        if m ==n
            Zmn(m,n) = 2*log(dx/a);
        else
%         p = (xb-xm)+sqrt((xb-xm).^2-a^2);   % Caluclation of solution matrix num
%         q = (xa-xm)+sqrt((xa-xm).^2-a^2);   % Caluclation of solution matrix den
%         Zmn(m,n) = log(abs(p/q));           % Getting Zmn matrix
%         A(i, j)=DELTA/abs(Y(i)-Y(j) )
        Zmn(m,n) = dx/(abs(xm(n)-xm(m)));
        end
    end
end
% Formation of Zmn matrix
% for m = 1:M
%     for n = 1:N
%         Zmn1(m,n) = Zmn(m,n);
%     end
% end
bm = 4*pi*ep_0*ones(M,1);           % Constant Matrix
an = Zmn\bm;
f1 = 0;
fvar = 0;
for n = 1:N
    xn = lim1+(n-1)*hx;
    if n == 1
        fvar = an(n)*rectpuls(x-xn,hx);
    else
    fvar = an(n)*rectpuls(x-xn,hx);
    end
    f1 = f1+fvar;
end
figure();
plot(x,f1*1e12,'b--','LineWidth', 2);
title('Densidade de Carga em relação a Seção do Fio');
xlabel('Seção do Fio');
ylabel('Densidade de Carga(10^-^1^0C/m^3)');
