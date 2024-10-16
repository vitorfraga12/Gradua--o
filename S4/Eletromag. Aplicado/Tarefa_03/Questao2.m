L=200;
U=zeros(L,L);
U(end/2-40,end/2-50:end/2+50)= 220 ;
U(end/2+40,end/2-50:end/2+50)=-220 ;
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
title (' Electric potential of the parralel plate Capacitor');
view(-54,6);
%  Electrostatic field
[Ex,Ey]=gradient(U);
figure, contour(U,'LineWidth',2);
hold on, quiver(Ex,Ey,4), hold off
title(' Electrical field of the parallel plate Capacitor in static case');
%Thrid figure for better viusalization of the electric field
figure,contour(U,'LineWidth',2); hold on, quiver(Ex,Ey,4), hold off;
zoom(3);
axis([10 40 130 160]);
xlabel(' Lx'), ylabel('Ly'), title(' Electrcial field :Zoom x4 Left Corner of the 1st plate ')
%% 
% Initial conditions
L=200;
U=zeros(L,L);
d = L/200; % 1m em escala da grade
U(:,1) = 220; % Fio na posição x=0 com potencial de 220V
U(:,1+d) = -220; % Fio deslocado em 1m com potencial de -220V

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
title (' Electric potential of the parallel wires');
view(-54,6);

% Electrostatic field
[Ex,Ey]=gradient(U);
figure, contour(U,'LineWidth',2);
hold on, quiver(Ex,Ey,4), hold off
title(' Electrical field of the parallel wires in static case');

% Third figure for better visualization of the electric field
figure, contour(U,'LineWidth',2); hold on, quiver(Ex,Ey,4), hold off;
zoom(3);
axis([10 40 130 160]);
xlabel(' Lx'), ylabel('Ly'), title(' Electrical field: Zoom x4 Near the wires ')
%% 
% Solving numerically the 2D Laplace Equation for parallel wires
% using finite differences method, convergence is attained using the norm's
% criterion with tolerance=6.00. Number of iteration N=611.
%
%  - Laplace eq : d²U(x,y)/dx²+d²U(x,y)/dy²=0
%  - boundaries : U(x=0,y)=0, U(x=L,0)=0, U(x,y=0)=0, U(x,y=L)=0.
%
% Parametrs     :
%  - Dimensions : square box of length L=200 mm .
%  - Voltage    : two wires : (1) at 220 volts and (2) at -220 volts.
%  - distance   : between wires d=1 m.
%  
%Initial conditions
L=200;
U=zeros(L,L);
d = 100; % Distance between the wires in grid points (1 meter = 100 grid points)
U(1:L,L/2) = 220; % Positive potential for left wire
U(1+d:L+d,L/2) = -220; % Negative potential for right wire

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
title (' Electric potential of the parallel wires');
view(-54,6);
%  Electrostatic field
[Ex,Ey]=gradient(U);
figure, contour(U,'LineWidth',2);
hold on, quiver(Ex,Ey,4), hold off
title(' Electrical field of the parallel wires in static case');
%Thrid figure for better viusalization of the electric field
figure,contour(U,'LineWidth',2); hold on, quiver(Ex,Ey,4), hold off;
zoom(3);
axis([10 40 130 160]);
xlabel(' Lx'), ylabel('Ly'), title(' Electrcial field :Zoom x4 Left side of the 1st wire ')

