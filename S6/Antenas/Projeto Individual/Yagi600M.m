%% Propriedades da Antena

% Criação da antena Yagi-Uda ajustada para 600 MHz
antennaObject = design(yagiUda, 600 * 1e6); % Define a frequência de operação em 600 MHz
antennaObject.Exciter = dipoleFolded; % Define o dipolo dobrado como elemento ativo

% Configuração das dimensões do dipolo dobrado
antennaObject.Exciter.Length = 0.21775; % Comprimento do dipolo em metros
antennaObject.Exciter.Width = 0.0055462; % Largura do dipolo em metros
antennaObject.Exciter.Spacing = 0.002848; % Espaçamento entre os condutores do dipolo dobrado

% Configuração do número de diretores
antennaObject.NumDirectors = 9; % Define o número de diretores como 9

% Visualizar a geometria da antena Yagi-Uda
figure;
show(antennaObject); % Mostra a estrutura física da antena
title('Geometria da Antena Yagi-Uda (600 MHz)');

%% Análise da Antena

% Definição da frequência para análise
plotFrequency = 600 * 1e6; % Frequência central para os gráficos (600 MHz)
freqRange = (540:6:660) * 1e6; % Faixa de frequência de 540 MHz a 660 MHz

% Impedância de referência do sistema
refImpedance = 50; % Impedância típica de sistemas de telecomunicação (50 ohms)

% Cálculo e plotagem da impedância em função da frequência
figure;
impedance(antennaObject, freqRange);
title('Impedância da Antena Yagi-Uda (600 MHz)');
xlabel('Frequência (Hz)');
ylabel('Impedância (Ohms)');

% Cálculo e plotagem dos parâmetros S (exemplo: S11)
figure;
s = sparameters(antennaObject, freqRange, refImpedance); % Calcula os parâmetros S
rfplot(s); % Plota o gráfico de S11
title('Parâmetro S11 da Antena Yagi-Uda (600 MHz)');
xlabel('Frequência (Hz)');
ylabel('|S11| (dB)');

% Padrão de radiação em 3D
figure;
pattern(antennaObject, plotFrequency); % Exibe o padrão de radiação tridimensional
title('Padrão de Radiação 3D da Antena Yagi-Uda (600 MHz)');

% Padrão de radiação no plano azimutal (horizontal)
figure;
patternAzimuth(antennaObject, plotFrequency, 0, 'Azimuth', 0:5:360); % Padrão de radiação no plano H
title('Padrão de Radiação no Plano Azimutal (Horizontal) - 600 MHz');
xlabel('Ângulo Azimutal (graus)');
ylabel('Ganho (dB)');

% Padrão de radiação no plano de elevação (vertical)
figure;
patternElevation(antennaObject, plotFrequency, 0, 'Elevation', 0:5:360); % Padrão de radiação no plano E
title('Padrão de Radiação no Plano de Elevação (Vertical) - 600 MHz');
xlabel('Ângulo de Elevação (graus)');
ylabel('Ganho (dB)');
