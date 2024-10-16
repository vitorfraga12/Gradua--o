a = 0; %Limite infeior
b = 2; %Limite superior 
h = (b-a)/2;
x1 = a;
x2 = a+h;
x3 = b;
f1 = 1/sqrt(4 + x1^2);
f2 = 1/sqrt(4 + x2^2);
f3 = 1/sqrt(4 + x3^2);
integral = (h/3)*(f1+ 4*f2 + f3);
%% 1/3 de Simpsons para Eletromag Aplicado
a = 0;
b = 5;
f = @(x) 1/sqrt(b^2 +x.^2);
n = 6;
h = (b-a)/n;
k = 1:1:n-1;
S = zeros(size(k)); % Inicialize S como um vetor de zeros

for i = 1:length(k)
    x = a + k(i) * h;
    S(i) = f(x);
end

Se = sum(S(2:2:end));
So = sum(S(1:2:end));
out = (h/3)*(f(a)+f(b)+(2.*So) + (4.*Se));