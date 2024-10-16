% %Projeto de filtro maior banda passante em guia retangular
% 
% %Parâmetros do filtro e do guia
% f = 5e9; % Frequência em Hz
% lb = 100e6; % Largura de banda em Hz
% 
% cg = 0.67; % Comprimento da guia em metros
% lg = 0.12; % Largura da guia em metros
% 
% %Variáveis
% eps = 8.854e-12;
% mu = 4*pi*1e-7;
% tt = 1e-9; % Tempo total em segundos
% dx = 0.01; % Passo de discretização de x em metros
% dy = 0.01; % Passo de discretização de y em metros
% dt = 0.25*dx/3e8; % Passo de discretização temporal em segundos
% 
% % Estruturação do guia de ondas
% x = -cg/2:dx:cg/2;
% y = -lg/2:dy:lg/2;
% [X,Y] = meshgrid(x,y);
% guia = zeros(size(X));
% guia(abs(X) <= cg/2 & abs(Y) <= lg/2) = 1;
% 
% Ez = zeros(size(X));
% Hy = zeros(size(X));
% ts = round(tt/dt);
% 
% for n = 1:ts
%     % Atualização dos campos elétricos (Ez)
%     Ez(:,2:end) = Ez(:,2:end) + (dt/(eps*dx))*(Hy(:,2:end) - Hy(:,1:end-1));
%     
%     % Atualização dos campos magnéticos (Hy)
%     Hy(:,1:end-1) = Hy(:,1:end-1) + (dt/(mu*dx))*(Ez(:,2:end) - Ez(:,1:end-1));
%     
%     % Aplicação da fonte Gaussiana
%     Ez((size(X,1)-1)/2,(size(X,2))/2) =exp(-((n-30)^2)/(100));
% end
% 
% subplot(2,1,1)
% % Plotagem da guia retangular
% imagesc(x, y, guia);
% colormap("colorcube");
% axis equal tight;
% xlabel('x em metros');
% ylabel('y em metros');
% title('Guia de Ondas Retangular');
% 
% % Plotagem do Campo Elétrico Ez
% subplot(2,1,2)
% imagesc(x, y, Ez);
% colormap("prism");
% axis equal tight;
% xlabel('x em metros');
% ylabel('y em metros');
% title('Campo Elétrico Ez na Guia Retangular');

% Projeto de filtro maior banda passante em guia circular

% Parâmetros do filtro e do guia
f = 5e9; % Frequência em Hz
lb = 100e6; % Largura de banda em Hz

rg = 0.4; % Raio da guia circular em metros
pnts = 500; % Número de pontos da guia circular

% Estruturação do guia circular
th = linspace(0, 2*pi, pnts);
x = rg*cos(th);
y = rg*sin(th);

%  Resposta em frequência do filtro para guia circular
fq = linspace(f - lb/2, f + lb/2, 1000);
rps = abs(sinc((fq - f) / lb));


subplot(1,2,1)
% Plotagem do guia circular
plot(x, y, 'k', 'LineWidth', 2);
axis equal;
title('Guia de Ondas Circular');
xlabel('X em metros');
ylabel('Y em metros');
grid on;

subplot(1,2,2)
% Plotagem da resposta em frequência
plot(fq, rps, 'k', 'LineWidth', 2);
title('Resposta em Frequência');
xlabel('Frequência em Hz');
ylabel('Magnitude');
grid on