import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erfc


def lin2dB(x):
    return 10*np.log10(x) # Converte a escala linear para dB

def dB2lin(x):
    return 10**(x/10) # Converte a escala dB para linear

def Qfunc(x):
    return 0.5*erfc(x/np.sqrt(2)) # Gera a função Q a partir da função erro complementar


####################################################################################################################################################################################

# Função para gerar a codificação de Gray
def getGrayCode(M_order):
    """
    Gera o código Gray para M níveis.
    
    Parâmetros:
    - M_order (int): Número de níveis na constelação.
    
    Retorna:
    - gray_codes (dict): Mapeamento de números binários para Gray.
    """
    if not (M_order & (M_order - 1)) == 0:  # Verifica se M é potência de 2
        raise ValueError("M deve ser uma potência de 2.")
    
    gray_codes = {} # Dicionário para mapear binário para Gray
    for i in range(M_order): 
        # Para fazer a conversão de binário para Gray, basta fazer um XOR entre o número e o número deslocado uma posição para a direita (>> 1) 
        gray = i ^ (i >> 1)  # Cálculo do código Gray
        gray_codes[i] = gray 
    
    return gray_codes

#####################################################################################################################################

# Função para gerar a constelação M-PAM
def getPAM(bits:np.ndarray, m_order: int, c_distance: float):
    ''' 
    Realiza a modulação M-PAM com codificação Gray.

    Parâmetros:
    - bits (array): Sequência de bits a ser modulada.
    - m_order (int): Número de níveis na constelação (M-PAM).
    - c_distance (float): Distância mínima entre níveis.

    Retorna:
    - symbols (array): Símbolos modulados em M-PAM. 
    '''
    if not np.log2(m_order).is_integer():
        raise ValueError('O número de níveis da constelação deve ser uma potência de 2.')
    
    # Definindo o número de bits por símbolo
    K_BitsPerSymbol = int(np.log2(m_order))
    if len(bits) % K_BitsPerSymbol != 0:
        raise ValueError("O número de bits deve ser múltiplo de log2(M).")
    
    # Gerando o código Gray
    gray_codes = getGrayCode(m_order)
    gray_map = {i: gray_codes[i] for i in range(m_order)}  # Índice para código Gray

    # Gerando os níveis da constelação (usando índices e escalando com c_distance)
    pam_levels = np.array([2 * i - (m_order - 1) for i in range(m_order)]) * c_distance

    # Mapeamento do código Gray para níveis da constelação
    mapping = {tuple(map(int, format(gray_map[i], f'0{K_BitsPerSymbol}b'))): pam_levels[i] for i in range(m_order)}

    # Modulação dos bits
    symbols = []
    for i in range(0, len(bits), K_BitsPerSymbol):
        bit_chunk = tuple(bits[i:i + K_BitsPerSymbol])
        if bit_chunk not in mapping:
            raise KeyError(f"Bit chunk {bit_chunk} não encontrado no mapeamento.")
        symbols.append(mapping[bit_chunk])


    return np.array(symbols)

#####################################################################################################################################

# Função para gerar a constelação M-QAM
def getQAM(bits: np.ndarray, m_order: int, c_distance: float):
    ''' 
    Realiza a modulação M-QAM com codificação Gray.

    Parâmetros:
    - bits (array): Sequência de bits a ser modulada.
    - m_order (int): Número de níveis na constelação (Q-PAM).
    - c_distance (float): Distância mínima entre níveis.

    Retorna:
    - symbols (array): Símbolos modulados em Q-PAM.
    - gray_labels (dict): Mapeamento de cada símbolo para o código Gray correspondente.
    '''
    if not np.log2(m_order).is_integer():
        raise ValueError("O número de níveis da constelação deve ser uma potência de 2.")
    
    sqrt_m = int(np.sqrt(m_order)) # Número de níveis em cada dimensão da constelação, que será utilizado para fazer √M-PAM
    if sqrt_m**2 != m_order:
        raise ValueError("O número de níveis da constelação deve ser um quadrado perfeito.")
    
    # Definindo o número de bits por símbolo e dimensão
    K_BitsPerDimention = int(np.log2(sqrt_m))
    K_BitsPerSymbol = 2*K_BitsPerDimention
    if len(bits) % K_BitsPerSymbol != 0:
        raise ValueError("O número de bits deve ser múltiplo de 2*log2(M).")
    
    symbols = []
    gray_labels = {}
    for i in range(0, len(bits), K_BitsPerSymbol):
        in_phase_bits = bits[i:i+K_BitsPerDimention]
        quadrature_bits = bits[i+K_BitsPerDimention:i+K_BitsPerSymbol]

        # Fazendo a modulação √M-PAM
        in_phase = getPAM(in_phase_bits, sqrt_m, c_distance)[0] # Modulação PAM na dimensão I (in-phase) da constelação 
        quadrature = getPAM(quadrature_bits, sqrt_m, c_distance)[0] # Modulação PAM na dimensão Q (quadrature) da constelação

        # Criando agora a constelação a representação complexa de cada símbolo
        symbol = in_phase + 1j*quadrature
        symbols.append(symbol)

        # Convertendo os bits para string e mapeando para o código Gray, a partir da concatenação dos bits de cada dimensão (I e Q)
        gray_labels[symbol] = ''.join(map(str, in_phase_bits)) + ''.join(map(str, quadrature_bits))

    return np.array(symbols), gray_labels

#####################################################################################################################################

# Função para a constelação M-PSK

def getPSK(bits: np.ndarray, m_order: int):
    """
    Realiza a modulação M-PSK com codificação Gray.

    Parâmetros:
    - bits (array): Sequência de bits a ser modulada.
    - m_order (int): Número de níveis na constelação (M-PSK).

    Retorna:
    - symbols (array): Símbolos modulados em M-PSK.
    - gray_labels (list): Lista de códigos Gray correspondentes a cada símbolo.
    """
    if not np.log2(m_order).is_integer():
        raise ValueError("O número de níveis da constelação deve ser uma potência de 2.")
    
    # Número de bits por símbolo
    K_BitsPerSymbol = int(np.log2(m_order))
    if len(bits) % K_BitsPerSymbol != 0:
        raise ValueError("O número de bits deve ser múltiplo de log2(M).")
    
    # Gerando o código Gray usando a função getGrayCode
    gray_codes = getGrayCode(m_order)
    
    # Gerando os símbolos M-PSK
    phases = np.linspace(0, 360, m_order, endpoint=False)  # Ângulos uniformemente distribuídos
    phases = np.radians(phases)  # Convertendo para radianos
    psk_symbols = np.exp(1j * phases)  # Símbolos no círculo unitário

    # Mapeando os códigos Gray para símbolos
    mapping = {gray_codes[i]: psk_symbols[i] for i in range(m_order)}

    # Modulação
    symbols = []
    gray_labels = {} 
    for i in range(0, len(bits), K_BitsPerSymbol):
        bit_chunk = int("".join(map(str, bits[i:i + K_BitsPerSymbol])), 2)  # Bits como inteiro
        gray_label = gray_codes[bit_chunk]  # Código Gray
        symbol = mapping[gray_label]  # Símbolo correspondente na constelação PSK
        
        symbols.append(symbol)
        gray_labels[symbol] = format(gray_label, f'0{K_BitsPerSymbol}b')  # Código Gray como string binária
    
    return np.array(symbols), gray_labels


#####################################################################################################################################

def getEnergy(symbols: np.ndarray):
    ''' 
    Calcula a energia média dos símbolos modulados.

    Parâmetros:
    - symbols (array): Símbolos modulados.

    Retorna:
    - energy (float): Energia média dos símbolos.
    '''
    
    return np.mean(np.abs(symbols)**2)

#####################################################################################################################################

def getMinDistance(symbols: np.ndarray):
    ''' 
    Calcula a menor distância entre os símbolos da constelação.

    Parâmetros:
    symbols (array): Símbolos modulados.

    Retorna:
    - min_distance (float): Menor distância entre os símbolos.
    '''
    # Vendo quais os símbolos únicos na constelação
    unique_symbols = np.unique(symbols)

    # Calcula a matriz de distâncias entre todos os pares de símbolos
    dist_matrix = np.abs(unique_symbols[:, None] - unique_symbols[None, :])
    
    # Exclui os valores diagonais (auto-distâncias iguais a 0)
    dist_matrix[dist_matrix == 0] = np.inf

    # Retorna a menor distância encontrada
    d_min = np.min(dist_matrix)
    return d_min

#####################################################################################################################################

def addAWGN(symbols: np.array, n_samples: int, snr: float, signal_energy: float):
    ''' 
    Adiciona ruído AWGN ao sinal modulado.

    Parâmetros:
    - symbols (array): Símbolos modulados.
    - n_samples (int): Número de amostras do ruído.
    - snr (float): Relação sinal-ruído em dB.
    - signal_energy (float): Energia média dos símbolos.

    Retorna:
    - tx_signal (array): Sinal modulado com ruído AWGN.
    '''
    # Calculando o desvio padrão do ruído
    noise_std_N0 = np.sqrt(signal_energy / (2 * dB2lin(snr))) # Desvio padrão do ruído
    noise = np.random.normal(0, noise_std_N0, n_samples) + 1j * np.random.normal(0, noise_std_N0, n_samples) # Gerando ruído
    tx_signal = symbols + noise # Adicionando o ruído ao sinal modulado
    return tx_signal

#####################################################################################################################################

def TheoricalErrorProbMQAM(M: int, snr: float):
    ''' 
    Calcula a probabilidade de erro teórica de símbolo para M-QAM.

    Parâmetros:
    - M (int): Número de níveis na constelação.
    - snr (float): Relação sinal-ruído em dB.

    Retorna:
    - error_prob (float): Probabilidade de erro de símbolo.
    '''
    sqrt_m = int(np.sqrt(M)) # Número de níveis em cada dimensão da constelação
    if sqrt_m**2 != M:
        raise ValueError("O número de níveis da constelação deve ser um quadrado perfeito.")
    
    # Cálculo da probabilidade de erro de símbolo para M-QAM
    snr_lin = dB2lin(snr)

    arg = np.sqrt((3*snr_lin)/(M-1))

    term_1 = 4 * ( 1 - (1/np.sqrt(M)) ) * Qfunc(arg)
    term_2 = 4 * ( 1 - (1/np.sqrt(M)) )**2 * Qfunc(arg)**2

    error_prob = term_1 - term_2
    return error_prob

#####################################################################################################################################

def TheoricalErrorProbMPSK(M: int, snr: float):
    ''' 
    Calcula a probabilidade de erro teórica de símbolo para M-PSK.

    Parâmetros:
    - M (int): Número de níveis na constelação.
    - snr (float): Relação sinal-ruído em dB.

    Retorna:
    - error_prob (float): Probabilidade de erro de símbolo.
    '''

    # Cálculo da probabilidade de erro de símbolo para M-PSK
    snr_lin = dB2lin(snr)

    arg = np.sqrt(2 * snr_lin) * np.sin(np.pi / M) # Argumento da Q-function

    error_prob = 2 * Qfunc(arg)

    return error_prob


#####################################################################################################################################

def slicer(rx_symbols: np.ndarray, constellation: np.ndarray):
    """
    Decodifica símbolos recebidos baseando-se na distância mínima até os símbolos da constelação.

    Parâmetros:
    - received_symbols (array): Símbolos recebidos no canal AWGN.
    - constellation (array): Símbolos da constelação M-QAM.

    Retorna:
    - decoded_symbols (array): Símbolos decodificados correspondentes aos mais próximos na constelação.
    """
    decoded_symbols = []
    for sym in rx_symbols:
        # Calcula a distância euclidiana entre o símbolo recebido e todos os símbolos da constelação
        distances = np.abs(sym - constellation)
        # Seleciona o símbolo da constelação com a menor distância
        min_index = np.argmin(distances)
        decoded_symbols.append(constellation[min_index])
    return np.array(decoded_symbols)


#####################################################################################################################################

def getSERandBER(bits: np.ndarray, M_order:int ,snr: np.ndarray, modulation_type: str ,c_distance=1):
    """
    Calcula a SER e BER simuladas para M-QAM em um canal AWGN.

    Parâmetros:
    - bits (array): Sequência de bits a ser transmitida.
    - M_order (int): Ordem da constelação M-QAM.
    - snr (array): Faixa de valores de SNR em dB.
    - c_distance (float): Distância mínima entre níveis da constelação.
    - modulation_type (str): Tipo de modulação (QAM, PSK, PAM).

    Retorna:
    - SER (list): Taxa de erro de símbolo simulada.
    - BER (list): Taxa de erro de bit simulada.
    """

    n_bits_per_symbol = int(np.log2(M_order))
    num_symbols = len(bits) // n_bits_per_symbol

    # Gerando a constelação com base no tipo de modulação
    if modulation_type.upper() == 'QAM':
        symbols, gray_labels = getQAM(bits, M_order, c_distance)
    elif modulation_type.upper() == 'PSK':
        symbols, gray_labels = getPSK(bits, M_order)
    else:
        raise ValueError("Tipo de modulação inválido. Use 'QAM' ou 'PSK'.")

    constellation = np.unique(symbols)
    signal_energy = getEnergy(symbols) # Calculando a energia média dos símbolos

    SER = []   
    BER = []

    for snr_value in snr:
        # Adicionando ruído AWGN ao sinal modulado
        symbols_rx = addAWGN(symbols, num_symbols, snr_value, signal_energy)

        # Decodificando os símbolos recebidos
        decoded_symbols = slicer(symbols_rx, constellation)

        # Obtendo a quantidade de erros de símbolo
        error_symbols = np.sum(symbols != decoded_symbols)
        SER.append(error_symbols / num_symbols)

        # Obtendo a quantidade de erros de bit
        error_bits = error_symbols * n_bits_per_symbol
        BER.append(error_bits / len(bits))

    return SER, BER


####################################################################################################################################################################################

def getSERandBER_defined(bits: np.ndarray, M_order:int ,snr: np.ndarray, constellation_values: np.ndarray):
    """
    Calcula a SER e BER simuladas para M-QAM em um canal AWGN.

    Parâmetros:
    - bits (array): Sequência de bits a ser transmitida.
    - M_order (int): Ordem da constelação M-QAM.
    - snr (array): Faixa de valores de SNR em dB.
    - c_distance (float): Distância mínima entre níveis da constelação.
    - constellation_values (array): Valores da constelação.

    Retorna:
    - SER (list): Taxa de erro de símbolo simulada.
    - BER (list): Taxa de erro de bit simulada.
    """

    n_bits_per_symbol = int(np.log2(M_order))
    num_symbols = len(bits) // n_bits_per_symbol

    symbols = constellation_values

    constellation = np.unique(symbols)
    signal_energy = getEnergy(symbols) # Calculando a energia média dos símbolos

    SER = []   
    BER = []

    for snr_value in snr:
        # Agora sim, passa os símbolos modulados corretamente para addAWGN
        symbols_rx = addAWGN(symbols, num_symbols, snr_value, signal_energy)


        # Decodificando os símbolos recebidos
        decoded_symbols = slicer(symbols_rx, constellation)

        # Obtendo a quantidade de erros de símbolo
        error_symbols = np.sum(symbols != decoded_symbols)
        SER.append(error_symbols / num_symbols)

        # Obtendo a quantidade de erros de bit
        error_bits = error_symbols * n_bits_per_symbol
        BER.append(error_bits / len(bits))

    return SER, BER, symbols_rx

####################################################################################################################################################################################

def map_region(x_i, x_j, prob_xi, prob_xj, var):
    """
    Detector de Máxima a Posteriori (MAP) para sinais recebidos.

    Parâmetros:
    - x_i (complex): Símbolo de comparação.
    - x_j (complex): Outro símbolo de comparação.
    - prob_xi (float): Probabilidade de ocorrência do símbolo x_i.
    - prob_xj (float): Probabilidade de ocorrência do símbolo x_j.
    - var (float): Variância do ruído.

    Retorna:
    - float: Valor da equação de decisão.
    """
    # Primeiro termo da equação obtida
    prim_term = (x_i**2 - x_j**2) / (2 * (x_i - x_j))
    
    # Segundo termo da equação obtida
    segun_term = (var / (x_i - x_j)) * (np.log(prob_xj / prob_xi))
    
    return prim_term + segun_term

####################################################################################################################################################################################

def apply_phase_offset(signal, phase_offset):
    '''
    Aplica um desvio de fase constante ao sinal recebido.

    Parâmetros:
    - signal (array): Sinal modulado recebido.
    - phase_offset (float): Desvio de fase em radianos.

    Retorna:
    - signal_with_phase (array): Sinal com desvio de fase aplicado.
    '''
    return signal * np.exp(1j * phase_offset)  # Multiplicação por um fator exponencial para rotacionar a fase

####################################################################################################################################################################################

def symbols_to_bits(symbols, constellation, m_order):
    '''
    Converte símbolos demodulados em uma sequência de bits utilizando o mapeamento Gray.

    Parâmetros:
    - symbols (array): Símbolos demodulados.
    - constellation (array): Conjunto de símbolos da constelação QPSK.
    - m_order (int): Ordem da modulação M-PSK.

    Retorna:
    - bits (array): Sequência de bits decodificados a partir dos símbolos.
    '''
    gray_codes = getGrayCode(m_order)  # Gera a sequência de códigos Gray
    bits_per_symbol = int(np.log2(m_order))  # Número de bits por símbolo
    bits = []

    for sym in symbols:
        distances = np.abs(sym - constellation)  # Calcula a distância euclidiana em relação à constelação
        closest_symbol_index = np.argmin(distances)  # Encontra o símbolo mais próximo
        gray_code = gray_codes[closest_symbol_index]  # Obtém o código Gray correspondente
        bits_chunk = list(format(gray_code, f'0{bits_per_symbol}b'))  # Converte o código Gray para uma sequência binária
        bits.extend([int(bit) for bit in bits_chunk])  # Adiciona os bits convertidos à lista

    return np.array(bits)  # Retorna a sequência de bits como um array

####################################################################################################################################################################################
def phase_correction(symbols_rx, constellation):
    """
    Estima e corrige o desvio de fase dos símbolos recebidos usando uma média ponderada para maior precisão.

    Parâmetros:
    - symbols_rx (array): Símbolos recebidos no canal.
    - constellation (array): Símbolos ideais da constelação.

    Retorna:
    - symbols_corrected (array): Símbolos corrigidos.
    - estimated_phase (float): Fase estimada no intervalo [0, 2π].
    """
    # Obtendo apenas os símbolos únicos da constelação
    unique_constellation = np.unique(constellation)

    # Obtendo as fases dos símbolos recebidos e da constelação
    phases_rx = np.angle(symbols_rx)
    phases_constellation = np.angle(unique_constellation)

    phase_differences = []
    weights = []

    # Comparação individual para evitar grandes matrizes (mais eficiente)
    for phase_rx in phases_rx:
        # Calcula a diferença de fase para cada símbolo recebido em relação à constelação
        phase_diff = np.abs(phase_rx - phases_constellation)
        phase_diff = (phase_diff + np.pi) % (2 * np.pi) - np.pi  # Ajusta para intervalo [-π, π]

        # Encontra a menor diferença de fase
        min_diff = np.min(phase_diff)
        phase_differences.append(min_diff)

        # Define um peso inversamente proporcional à diferença de fase (quanto menor a diferença, maior o peso)
        weight = np.exp(-min_diff)  # Ponderação exponencial para dar mais peso às estimativas confiáveis
        weights.append(weight)

    # Convertendo para numpy arrays
    phase_differences = np.array(phase_differences)
    weights = np.array(weights)

    # Normalizando os pesos
    weights /= np.sum(weights)

    # Calculando a estimativa da fase como uma média ponderada
    estimated_phase = np.sum(weights * phase_differences)

    # Ajusta para o intervalo [0, 2π]
    estimated_phase = (estimated_phase + 2 * np.pi) % (2 * np.pi)

    # Corrige a fase dos símbolos recebidos
    symbols_corrected = symbols_rx * np.exp(-1j * estimated_phase)

    return symbols_corrected, estimated_phase




########################################################################################################################################################

def getSERandBER_without_awgn(bits: np.ndarray, M_order:int ,snr: np.ndarray, constellation_values: np.ndarray):
    """
    Calcula a SER e BER simuladas para M-QAM em um canal sem ruído.

    Parâmetros:
    - bits (array): Sequência de bits a ser transmitida.
    - M_order (int): Ordem da constelação M-QAM.
    - snr (array): Faixa de valores de SNR em dB.
    - c_distance (float): Distância mínima entre níveis da constelação.
    - constellation_values (array): Valores da constelação.

    Retorna:
    - SER (list): Taxa de erro de símbolo simulada.
    - BER (list): Taxa de erro de bit simulada.
    """

    n_bits_per_symbol = int(np.log2(M_order))
    num_symbols = len(bits) // n_bits_per_symbol

    symbols = constellation_values

    constellation = np.unique(symbols)

    SER = []   
    BER = []

    for snr_value in snr:
        # Agora sim, passa os símbolos modulados corretamente para addAWGN
        symbols_rx = symbols

        # Decodificando os símbolos recebidos
        decoded_symbols = slicer(symbols_rx, constellation)

        # Obtendo a quantidade de erros de símbolo
        error_symbols = np.sum(symbols != decoded_symbols)
        SER.append(error_symbols / num_symbols)

        # Obtendo a quantidade de erros de bit
        error_bits = error_symbols * n_bits_per_symbol
        BER.append(error_bits / len(bits))

    return SER, BER, symbols_rx

####################################################################################################################################################################################

def getSERandBER_withphase(bits: np.ndarray, M_order:int ,snr: np.ndarray, modulation_type: str , phase: float,c_distance=1):
    """
    Calcula a SER e BER simuladas para M-QAM em um canal AWGN.

    Parâmetros:
    - bits (array): Sequência de bits a ser transmitida.
    - M_order (int): Ordem da constelação M-QAM.
    - snr (array): Faixa de valores de SNR em dB.
    - c_distance (float): Distância mínima entre níveis da constelação.
    - modulation_type (str): Tipo de modulação (QAM, PSK, PAM).
    - phase (float): Desvio de fase em radianos

    Retorna:
    - SER (list): Taxa de erro de símbolo simulada.
    - BER (list): Taxa de erro de bit simulada.
    """

    n_bits_per_symbol = int(np.log2(M_order))
    num_symbols = len(bits) // n_bits_per_symbol

    # Gerando a constelação com base no tipo de modulação
    if modulation_type.upper() == 'QAM':
        symbols, gray_labels = getQAM(bits, M_order, c_distance)
    elif modulation_type.upper() == 'PSK':
        symbols, gray_labels = getPSK(bits, M_order)
    else:
        raise ValueError("Tipo de modulação inválido. Use 'QAM' ou 'PSK'.")

    constellation = np.unique(symbols)
    signal_energy = getEnergy(symbols) # Calculando a energia média dos símbolos

    SER = []   
    BER = []

    for snr_value in snr:
        # Adicionando ruído AWGN ao sinal modulado
        symbols_rx = (addAWGN(symbols, num_symbols, snr_value, signal_energy))*np.exp(1j*phase)

        # Decodificando os símbolos recebidos
        decoded_symbols = slicer(symbols_rx, constellation)

        # Obtendo a quantidade de erros de símbolo
        error_symbols = np.sum(symbols != decoded_symbols)
        SER.append(error_symbols / num_symbols)

        # Obtendo a quantidade de erros de bit
        error_bits = error_symbols * n_bits_per_symbol
        BER.append(error_bits / len(bits))

    return SER, BER


def getSERandBER_withphase_estimation(bits: np.ndarray, M_order:int, snr: np.ndarray, modulation_type: str, phase: float, c_distance=1):
    """
    Calcula SER e BER aplicando correção de fase com o estimador definido.

    Parâmetros:
    - bits (array): Sequência de bits a ser transmitida.
    - M_order (int): Ordem da constelação M-PSK.
    - snr (array): Faixa de valores de SNR em dB.
    - modulation_type (str): Tipo de modulação (QAM, PSK, PAM).
    - phase (float): Desvio de fase em radianos.

    Retorna:
    - SER (list): Taxa de erro de símbolo simulada.
    - BER (list): Taxa de erro de bit simulada.
    """

    n_bits_per_symbol = int(np.log2(M_order))
    num_symbols = len(bits) // n_bits_per_symbol

    # Gerando a constelação com base no tipo de modulação
    if modulation_type.upper() == 'QAM':
        symbols, gray_labels = getQAM(bits, M_order, c_distance)
    elif modulation_type.upper() == 'PSK':
        symbols, gray_labels = getPSK(bits, M_order)
    else:
        raise ValueError("Tipo de modulação inválido. Use 'QAM' ou 'PSK'.")

    constellation = np.unique(symbols)
    signal_energy = getEnergy(symbols)

    SER = []
    BER = []

    for snr_value in snr:
        # Adicionando ruído AWGN com desvio de fase
        symbols_rx = addAWGN(symbols, num_symbols, snr_value, signal_energy) * np.exp(1j * phase)

        symbols_corrected, estimated_phase = phase_correction(symbols_rx, symbols)


        # Decodificando os símbolos corrigidos
        decoded_symbols = slicer(symbols_corrected, constellation)

        # Calculando SER
        error_symbols = np.sum(symbols != decoded_symbols)
        SER.append(error_symbols / num_symbols)

        # Calculando BER
        decoded_bits = symbols_to_bits(decoded_symbols, constellation, M_order)
        error_bits = np.sum(bits != decoded_bits)
        BER.append(error_bits / len(bits))

    return SER, BER
