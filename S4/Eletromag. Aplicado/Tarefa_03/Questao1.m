% Define constants
ER = 1.0;
EO = 8.8541e-12;
AA = 1.0;
BB = 1.0;
D = 1.0;
VO = 2.0;

% Vary the value of M from 1 to a max value, say 50
maxM = 50;
C_values = zeros(1, maxM);

for M = 1:maxM
    N = M^2;
    NT = 2*N;
    DX = AA/M;
    DY = BB/M;
    DL = DX;

    % Calculate the elements of the coefficient matrix A
    K = 0;
    X = zeros(1, NT);
    Y = zeros(1, NT);
    for K1 = 1:2
        for K2 = 1:M
            for K3 = 1:M
                K = K + 1;
                X(K) = DX*(K2 - 0.5);
                Y(K) = DY*(K3 - 0.5);
            end
        end
    end
    Z = [zeros(1, N), D * ones(1, N)];

    A = zeros(NT, NT);
    for I = 1:NT
        for J = 1:NT
            if(I == J)
                A(I,J) = DL*0.8814/(pi*EO);
            else
                R = sqrt((X(I)-X(J))^2 + (Y(I)-Y(J))^2 + (Z(I)-Z(J))^2);
                A(I,J) = DL^2/(4.*pi*EO*R);
            end
        end
    end

    % Determine the matrix of constant vector B
    B = [ones(1, N), -1 * ones(1, N)];

    % Invert A and calculate RHO
    F = inv(A);
    RHO = F*B';
    SUM = sum(RHO(1:N));

    Q = SUM*(DL^2);
    C_values(M) = abs(Q)/VO;
end

% Plot the capacitance against M
figure;
plot(1:maxM, C_values);
xlabel('N (Discretização)');
ylabel('Capacitância (F)');
title('Capacitância x Discretização');
grid on;
%% 
Nx = 101;     % Number of X-grids
Ny = 101;     % Number of Y-grids
mpx = ceil(Nx/2); % Mid-point of x
mpy = ceil(Ny/2); % Mid point of y
   
Ni = 750;  % Number of iterations for the Poisson solver
V = zeros(Nx,Ny);   % Potential (Voltage) matrix
T = 0;            % Top-wall potential
B = 0;            % Bottom-wall potential
L = 0;            % Left-wall potential
R = 0;            % Right-wall potential
%-------------------------------------------------------------------------%
% Initializing edges potentials
%-------------------------------------------------------------------------%
V(1,:) = L;
V(Nx,:) = R;
V(:,1) = B;
V(:,Ny) = T;
%-------------------------------------------------------------------------%
% Initializing Corner potentials
%-------------------------------------------------------------------------%
V(1,1) = 0.5*(V(1,2)+V(2,1));
V(Nx,1) = 0.5*(V(Nx-1,1)+V(Nx,2));
V(1,Ny) = 0.5*(V(1,Ny-1)+V(2,Ny));
V(Nx,Ny) = 0.5*(V(Nx,Ny-1)+V(Nx-1,Ny));
%-------------------------------------------------------------------------%
length_plate = 51;  % Length of plate in terms of number of grids  
lp = floor(length_plate/2);
position_plate = 15; % Position of plate on x axis
pp1 = mpx+position_plate;
pp2 = mpx-position_plate;
for z = 1:Ni    % Number of iterations
        
        for i=2:Nx-1
        for j=2:Ny-1      
            
            % The next two lines are meant to force the matrix to hold the 
            % potential values for all iterations
            
                V(pp1,mpy-lp:mpy+lp) = 100;
                V(pp2,mpy-lp:mpy+lp) = -100;
                
                V(i,j)=0.25*(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1));
        end
        end
        
end
% Take transpose for proper x-y orientation
V = V';
[Ex,Ey]=gradient(V);
Ex = -Ex;
Ey = -Ey;
% Electric field Magnitude
 E = sqrt(Ex.^2+Ey.^2);  
x = (1:Nx)-mpx;
y = (1:Ny)-mpy;
% Contour Display for electric potential
figure(1)
contour_range_V = -101:0.5:101;
contour(x,y,V,contour_range_V,'linewidth',0.5);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x-axis in meters','fontsize',14);
ylabel('z-axis in meters','fontsize',14);
title('Electric Potential distribution, V(x,z) in volts','fontsize',14);
h1=gca;
set(h1,'fontsize',14);
fh1 = figure(1); 
set(fh1, 'color', 'white')
% Contour Display for electric field
figure(2)
contour_range_E = -20:0.05:20;
contour(x,y,E,contour_range_E,'linewidth',0.5);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x (m)','fontsize',14);
ylabel('z (m)','fontsize',14);
title('Campo Elétrico no plano y = 0','fontsize',14);
h2=gca;
set(h2,'fontsize',14);
fh2 = figure(2); 
set(fh2, 'color', 'white')
% Quiver Display for electric field Lines
figure(3)
contour(x,y,E,'linewidth',0.5);
hold on, quiver(x,y,Ex,Ey,2)
title('Campo Elétrico no plano y = 0','fontsize',14);
axis([min(x) max(x) min(y) max(y)]);
colorbar('location','eastoutside','fontsize',14);
xlabel('x (m)','fontsize',14);
ylabel('z (m)','fontsize',14);
h3=gca;
set(h3,'fontsize',14);
fh3 = figure(3); 
set(fh3, 'color', 'white')
%-------------------------------------------------------------------------%
% REFERENCE
%           SADIKU, ELEMENTS OF ELECTROMAGNETICS, 4TH EDITION, OXFORD
%-------------------------------------------------------------------------%
%% 
L=200;
U=zeros(L,L);
U(end/2-40,end/2-80:end/2+80)= 220 ;
U(end/2+40,end/2-80:end/2+80)=-220 ;
count=1;
tolerance=6.00;
norm_diff=100*rand(1);
% loops
while norm_diff>tolerance
    Norm1=norm(U);
    for i=2:L-1
        for j=2:L-1
            % let the voltage's conditions intact
            if U(i,j)==220 || U(i,j)==-220
                continue;
            end
            % Laplac equation  : Delta U=0 the algorithm :
            U(i,j)=(U(i-1,j)+U(i+1,j)+U(i,j-1)+U(i,j+1))./4;       
        end
    end
    Norm2=norm(U);
    norm_diff=abs(Norm2-Norm1);
    count=count+1;
end
fprintf(' Number of iterations N=%d\n',count);
% plot the results
surf(U);shading interp; colorbar;
xlabel(' Lx');
ylabel(' Ly');
zlabel(' Potential, in Volts');
title ('Potencial (V) no plano z = 0');
view(-54,6);
%  Electrostatic field
[Ex,Ey]=gradient(U);
figure, contour(U,'LineWidth',2);
hold on, quiver(Ex,Ey,4), hold off
title('Campo Elétrico no plano z = 0');
%Thrid figure for better viusalization of the electric field
figure,contour(U,'LineWidth',2); hold on, quiver(Ex,Ey,4), hold off;
zoom(3);
axis([10 40 130 160]);
xlabel(' Lx'), ylabel('Ly'), title(' Electrcial field :Zoom x4 Left Corner of the 1st plate ')
