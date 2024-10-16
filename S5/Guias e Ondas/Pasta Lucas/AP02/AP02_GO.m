%QUESTÃO 01

%Corrente de Entrada da Linha (Is);

%Define os valores fixos que serão usados na questão;
V = 100;
f1 = 8e6;
f2 = 10e6;
f = f1:2:f2;
Z0 = 40.82;
d = 32;
bt = 1.54e-8*f;

%Coeficiente de Reflexão da Fonte;
CRF = -0.0768*exp(-j*2*d.*bt);

%Corrente de Entrada da Linha;
CEL = (V.*(1-CRF))./(Z0.*(1+CRF));

%Corrente da Carga;
CC = (V.*(exp(-j*d.*bt)).*(1.0768))./(Z0.*(1 + CRF));

%Resposta em Magnitude e Fase para a Corrente de Entrada da Linha;
magCEL = abs(CEL);
argCEL = angle(CEL);

%Resposta em Magnitude e Fase da Corrente da Carga;
magCC = abs(CC);
argCC = angle(CC);

%Plot das Respostas em Magnitude e Fase para a Corrente de Entrada de Linha;
figure
subplot(2,1,1)
plot(f,magCEL,LineWidth=2,Color=[0 0 0])
title("Resposta em Frequência da Magnitude de Corrente de Entrada da Linha")
xlabel("Frequência (Hz)")
ylabel("Magnitude de CEL")
grid on
subplot(2,1,2)
plot(f,argCEL,LineWidth=2,Color=[0 0 0])
title("Resposta em Frequência da Fase de Corrente de Entrada da Linha")
xlabel("Frequência (Hz)")
ylabel("Fase de CEL")
grid on

%Plot das Respostas em Magnitude e Fase para a Corrente da Carga;
figure
subplot(2,1,1)
plot(f,magCC,LineWidth=2,Color=[0 0 0])
title("Resposta em Frequência da Magnitude de Corrente da Carga")
xlabel("Frequência (Hz)")
ylabel("Magnitude de CC")
grid on
subplot(2,1,2)
plot(f,argCC,LineWidth=2,Color=[0 0 0])
title("Resposta em Frequência da Fase de Corrente da Carga")
xlabel("Frequência (Hz)")
ylabel("Fase de CC")
grid on
