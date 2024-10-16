%% Questão 05
% Carregamos e armazenamos a imagem
imagem = imread('Lenna.png');

% Dimensões da imagem
[Nv, Nh, ~] = size(imagem);

% Sinais unidimensionais xh (horizontal) e xv (vertical)
nV = 1:Nv; % Amostras na direção vertical
nH = 1:Nh; % Amostras na direção horizontal

% a) x[-nV, nH]
xh_a = imagem(Nv:-1:1, nH, :); % Começamos nV pelo final para representar -nV

% b) x[nV, -nH]
xv_b = imagem(nV, Nh:-1:1, :); % Começamos nH pelo final para representar -nH

% c) x[-nV, -nH]
xh_c = imagem(Nv:-1:1, Nh:-1:1, :); % Começamos nH e nV pelo final para representar -nH e -nV

% d) x[nV - n0, nH], onde n0 é um valor inteiro
n0 = 50; % Valor que pode ser alterado
xh_d = imagem(nV(1:end-n0), nH, :);

% e) x[nV, nH - n1], onde n1 é um valor inteiro
n1 = 30; % Valor que pode ser alterado
xv_e = imagem(nV, nH(1:end-n1), :);

% f) x[nV - n2, nH - n3], onde n2 e n3 são valores inteiros, podemos
% substituir também por n0 e n1 dos itens anteriores
n2 = 20; % Valor que pode ser alterado
n3 = 10; % Valor que pode ser alterado
xh_f = imagem(nV(1:end-n2), nH(1:end-n3), :);

% Plotamos a imagem
figure;

subplot(3, 3, 1);
imshow(imagem);
title('Imagem Original');

subplot(3, 3, 2);
imshow(uint8(xh_a));
title('a) x[-nV, nH]');

subplot(3, 3, 3);
imshow(uint8(xv_b));
title('b) x[nV, -nH]');

subplot(3, 3, 4);
imshow(uint8(xh_c));
title('c) x[-nV, -nH]');

subplot(3, 3, 5);
imshow(uint8(xh_d));
title(['d) x[nV - ' num2str(n0) ', nH]']);

subplot(3, 3, 6);
imshow(uint8(xv_e));
title(['e) x[nV, nH - ' num2str(n1) ']']);

subplot(3, 3, 7);
imshow(uint8(xh_f));
title(['f) x[nV - ' num2str(n2) ', nH - ' num2str(n3) ']']);