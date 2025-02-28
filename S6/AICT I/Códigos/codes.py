import numpy as np
import random
from scipy.signal import find_peaks

def lin2db(x):
    return 10*np.log10(x)

def db2lin(x):
    return 10**(x/10)

def gerate_a_ula(m_antennas:int, d_in:float, theta_i:float):
    '''
    Função que gera a matriz A_{ula} 
    m_antennas (int): número de antenas do arranjo
    d_in (float): número de direções de chegada
    theta_i (float): ângulo de chegada em graus

    return: A_{ula} (np.array): matriz A_{ula}
    
    '''



    mu_spatial_frequency = -np.pi*np.sin(np.radians(theta_i))

    A_ula = np.zeros((m_antennas, d_in), dtype=complex)

    for col in range (d_in):
        for row in range (m_antennas):
            A_ula[row, col] = np.exp(1j * row * mu_spatial_frequency[col])

    
    return A_ula


def generate_signalawgn(A_ula:np.ndarray ,m_antennas:int, arrival_distance: int, t_snapshot: int, snr:float):
    """
    Gera um sinal recebido em um arranjo de antenas com ruído AWGN.
    
    Parâmetros:
    - A_ula (np.ndarray): Matriz de resposta direcional.
    - m_antennas (int): Número de antenas no arranjo.
    - arrival_distance (int): Número de sinais incidentes (direções de chegada).
    - t_snapshot (int): Número de snapshots.
    - snr (float): Relação sinal-ruído (em dB).
    
    """

    # Gerar o sinal transmitido
    sinal = np.zeros((arrival_distance, t_snapshot), dtype=complex)
    for i in range(arrival_distance):
        sinal[i] =((np.random.normal(size=t_snapshot) + 1j * np.random.normal(size=t_snapshot))/arrival_distance) / np.sqrt(2)

    # Sinal recebido no arranjo de antenas
    sinal_aula = np.dot(A_ula, sinal)

    # Adicionando ruído branco gaussiano ao sinal
    ruido_amplitude = 1/np.sqrt(db2lin(snr))  # Ajusta o ruído com base na SNR
    noise = (np.random.normal(0, ruido_amplitude, (m_antennas, t_snapshot)) + 
             1j * np.random.normal(0, ruido_amplitude, (m_antennas, t_snapshot))) / np.sqrt(2)
    sinal_final = sinal_aula + noise

    return sinal_final



def generate_signal(m_antennas:int, noise_subspace):
    '''
    Função que gera o espectro de potência do sinal recebido
    m_antennas (int): número de antenas do arranjo
    noise_subspace (np.array): subespaço de ruído

    return: angles (np.array): ângulos avaliados, p_spectrum (np.array): espectro de potência do sinal recebido.
    
    '''
    
    # Estimando o espectro de potência
    angles = np.linspace(-90, 90, 181)
    p_spectrum = np.zeros(angles.shape)

    # Calculando o espectro de potência
    for index_angle, angle in enumerate(angles):
        steering_vector = np.exp(1j * np.arange(m_antennas) * (-np.pi * np.sin(np.radians(angle)))) # Vetor de direção, que é o vetor de entrada da matriz A

        numerator = np.abs(np.dot(np.conj(steering_vector.T), steering_vector))
        denominator = np.abs(np.dot(np.conj(steering_vector.T), np.dot(noise_subspace, np.dot(np.conj(noise_subspace.T), steering_vector))))
        
        p_spectrum[index_angle] = numerator / denominator if denominator != 0 else 0

    return angles, p_spectrum



def generate_music(A_ula: np.ndarray, arrival_distance: int, t_snapshot: int, m_antennas: int, snr: float):
    """
    Executa o algoritmo MUSIC para análise espectral.
    
    Parâmetros:
    - A_ula (np.ndarray): Matriz de resposta direcional.
    - arrival_distance (int): Número de sinais incidentes (direções de chegada).
    - t_snapshot (int): Número de snapshots.
    - m_antennas (int): Número de antenas no arranjo.
    - snr (float): Relação sinal-ruído (em dB).
    
    Retorno:
    - angles (np.ndarray): Ângulos (em graus) avaliados.
    - p_spectrum (np.ndarray): Espectro MUSIC correspondente aos ângulos.
    """

    # Gerar o sinal recebido com ruído AWGN
    sinal_final = generate_signalawgn(A_ula, m_antennas, arrival_distance, t_snapshot, snr)

    # Matriz de autocorrelação
    sinal_final_hermetiano = np.conj(sinal_final.T)
    autocor_matrix_estimada = np.dot(sinal_final, sinal_final_hermetiano) / t_snapshot

    # Decomposição em autovalores e autovetores
    eigenvalues, eigenvec = np.linalg.eig(autocor_matrix_estimada)

    # Ordenando os autovalores em ordem decrescente
    idx = eigenvalues.argsort()[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvec = eigenvec[:, idx]

    # Subespaço de ruído
    noise_subspace = eigenvec[:, arrival_distance:]

    # Gerar o espectro MUSIC
    angles, p_spectrum = generate_signal(m_antennas, noise_subspace)

    return angles, p_spectrum



def generate_angles_with_min_diff(d_arrival, min_diff, lower=-50, upper=50):
    """
    Gera ângulos uniformemente distribuídos com uma diferença mínima entre eles.
    
    Parâmetros:
    - d_arrival (int): Número de ângulos a serem gerados.
    - min_diff (float): Diferença mínima permitida entre ângulos consecutivos.
    - lower (float): Limite inferior do intervalo de geração.
    - upper (float): Limite superior do intervalo de geração.
    
    Retorno:
    - np.ndarray: Ângulos gerados com diferença mínima garantida. (em graus)
    """
    while True:
        phi_uniform = np.random.uniform(lower, upper, d_arrival)
        phi_uniform_sorted = np.sort(phi_uniform)
        if np.all(np.diff(phi_uniform_sorted) >= min_diff):
            return phi_uniform_sorted


def find_peaks_d(p_spectrum, d_arrival, prominence=2): 
    '''
    Encontra os picos no espectro MUSIC, considerando a proeminência dos picos.
    
    Parâmetros:
    - p_spectrum (np.ndarray): Espectro MUSIC.
    - d_arrival (int): Número de picos desejados.
    - prominence (float): Proeminência mínima para considerar um pico.
    
    Retorno:
    - peaks (np.ndarray): Índices dos picos encontrados.
    - values (np.ndarray): Valores dos picos encontrados.
    '''
    # Encontre os picos com base na proeminência
    peaks, properties = find_peaks(p_spectrum, prominence=prominence)

    # Classificar os picos pela proeminência (valores mais altos)
    locs = np.argsort(properties["prominences"])[::-1][:d_arrival]  # Ordena pelos maiores picos de proeminência
    locs = peaks[locs]  # Obtém os índices dos picos mais proeminentes
    locs = locs - 90
    locs = np.sort(locs)
    values = p_spectrum[locs]  # Valores dos picos mais proeminentes

    return locs, values


def find_rmse(snr_values: np.ndarray, iterations: int, m_antennas: int, d_arrival: int, t_snapshot: int):
    '''
    Calcula o RMSE (Root Mean Square Error) dos ângulos estimados pelo algoritmo MUSIC em diferentes valores de SNR.

    Parâmetros:
    - snr_values (list or np.ndarray): Lista de valores de SNR (em dB) a serem avaliados.
    - iterations (int): Número de iterações por valor de SNR.
    - m_antennas (int): Número de antenas no arranjo.
    - d_arrival (int): Número de ângulos de chegada (direções de chegada).
    - t_snapshot (int): Número de snapshots para estimativa de autocorrelação.

    Retorno:
    - rmse_maior (list): Lista com os valores de RMSE para os maiores ângulos estimados, correspondentes a cada SNR.
    - rmse_menor (list): Lista com os valores de RMSE para os menores ângulos estimados, correspondentes a cada SNR.
    '''
    
    rmse_maior = []
    rmse_menor = []
    
    for snr_index in snr_values:
        phi_maior = []
        phi_menor = []

        for iteration_index in range(iterations):
            # Gerando os ângulos de chegada
            phit_uniform = generate_angles_with_min_diff(d_arrival, 20)

            A_ula = gerate_a_ula(m_antennas, d_arrival, phit_uniform)
            angles, p_spectrum = generate_music(A_ula, d_arrival, t_snapshot, m_antennas, snr_index)

            # Calculando os ângulos de pico estimados pelo MUSIC
            top_peak_angles, _ = find_peaks_d(p_spectrum, d_arrival)

    

            maior_diferenca = phit_uniform[0] - top_peak_angles[0] 
            menor_diferena = phit_uniform[1] - top_peak_angles[1]

            phi_maior.append(maior_diferenca)
            phi_menor.append(menor_diferena)

        # Calculando o RMSE para os ângulos estimados
        phi_maior = np.array(phi_maior)
        phi_menor = np.array(phi_menor)

        rmse_maior.append(np.sqrt(np.mean(np.abs(phi_maior) ** 2)))
        rmse_menor.append(np.sqrt(np.mean(np.abs(phi_menor) ** 2)))

    return rmse_maior, rmse_menor

def root_music(signal_awgn: np.ndarray, m_antennas: int, d_arrival: int, true_angles: np.ndarray):
    """
    Implementação do algoritmo Root-MUSIC para estimação de ângulos de chegada (DoA).

    Parâmetros:
    - signal_awgn (np.ndarray): Sinais recebidos com ruído.
    - m_antennas (int): Número de antenas.
    - d_arrival (int): Número de ângulos de chegada (fontes).
    - true_angles (np.ndarray): Ângulos verdadeiros para comparação.

    Retorno:
    - estimated_angles (np.ndarray): Ângulos estimados pelo Root-MUSIC.
    - selected_roots (np.ndarray): Raízes selecionadas para conversão em ângulos.
    """

    # Parâmetros da onda
    wavelength = 1  # Comprimento de onda normalizado
    d = wavelength / 2  # Espaçamento entre sensores
    k = 2 * np.pi / wavelength  # Número de onda

    # Decomposição SVD
    U, S, Vh = np.linalg.svd(signal_awgn)
    noise_subspace = U[:, d_arrival:]  # Subespaço de ruído

    # Construção da matriz C, definida como C = E_n E_n^H, onde E é o subespaço do ruído
    C = np.dot(noise_subspace,noise_subspace.conj().T)  

    #  Definimos a ordem do polinômio como (2M-1), sendo M o número de elementos de antenas
    z = np.zeros(2 * m_antennas - 1, dtype=complex)

    # Soma das diagonais da matriz C
    for l in range(m_antennas):
        z[l] = np.sum(np.diag(C, l))
        z[-(l + 1)] = z[l].conj()  # Mantendo a simetria

    # Encontrando as raízes do polinômio
    z_roots = np.roots(z)

    # Selecionando as raízes mais próximas do círculo unitário
    z_roots = z_roots[np.abs(np.abs(z_roots) - 1) < 0.1]  

    # Convertendo raízes para ângulos
    angles = np.arcsin(np.angle(z_roots) / (k * d))*180/np.pi
    
    # Seleção das raízes mais próximas dos ângulos reais

    selected_roots = []
    selected_angles = []
    for true_angle in true_angles:
        idx = np.argmin(np.abs(angles - true_angle))
        selected_roots.append(z_roots[idx])
        selected_angles.append(angles[idx])

    return np.array(selected_angles), np.array(selected_roots)

def find_rmse_rootmusic(snr_values: np.ndarray, iterations: int, m_antennas: int, d_arrival: int, t_snapshot: int):
    '''
    Calcula o RMSE (Root Mean Square Error) dos ângulos estimados pelo algoritmo Root MUSIC em diferentes valores de SNR.

    Parâmetros:
    - snr_values (list or np.ndarray): Lista de valores de SNR (em dB) a serem avaliados.
    - iterations (int): Número de iterações por valor de SNR.
    - m_antennas (int): Número de antenas no arranjo.
    - d_arrival (int): Número de ângulos de chegada (direções de chegada).
    - t_snapshot (int): Número de snapshots para estimativa de autocorrelação.

    Retorno:
    - rmse_maior (list): Lista com os valores de RMSE para os maiores ângulos estimados, correspondentes a cada SNR.
    - rmse_menor (list): Lista com os valores de RMSE para os menores ângulos estimados, correspondentes a cada SNR.
    '''
    
    rmse_maior = []
    rmse_menor = []
    
    for snr_index in snr_values:
        phi_maior = []
        phi_menor = []

        for iteration_index in range(iterations):
            # Gerando os ângulos de chegada
            phit_uniform = generate_angles_with_min_diff(d_arrival, 20)

            A_ula = gerate_a_ula(m_antennas, d_arrival, phit_uniform)

            # Geração dos sinais recebidos
            signal_awgn = generate_signalawgn(A_ula, m_antennas, d_arrival, t_snapshot, snr_index)

            # Cálculo do Root-MUSIC
            estimated_angles, selected_roots = root_music(signal_awgn, m_antennas, d_arrival, phit_uniform)

            maior_diferenca = phit_uniform[0] - estimated_angles[0] 
            menor_diferena = phit_uniform[1] -  estimated_angles[1]

            phi_maior.append(maior_diferenca)
            phi_menor.append(menor_diferena)

        # Calculando o RMSE para os ângulos estimados
        phi_maior = np.array(phi_maior)
        phi_menor = np.array(phi_menor)

        rmse_maior.append(np.sqrt(np.mean(np.abs(phi_maior) ** 2)))
        rmse_menor.append(np.sqrt(np.mean(np.abs(phi_menor) ** 2)))

    return rmse_maior, rmse_menor

def mvdr_beamforming(A_ula: np.ndarray, m_antennas: int, arrival_distance: int, t_snapshot: int, snr: float, theta_target: float):
    """
    Implementação do beamforming MVDR.

    Parâmetros:
    - A_ula (np.ndarray): Matriz de resposta direcional.
    - m_antennas (int): Número de antenas no arranjo.
    - arrival_distance (int): Número de ângulos de chegada (fontes).
    - t_snapshot (int): Número de snapshots.
    - snr (float): Relação sinal-ruído (em dB).
    - theta_target (float): Ângulo de interesse (graus).

    Retorno:
    - B (np.ndarray): Padrão de radiação (beampattern).
    - theta_scan (np.ndarray): Ângulos varridos (-90° a 90°).
    """

    # Geração dos sinais recebidos com ruído
    received_signal = generate_signalawgn(A_ula, m_antennas, arrival_distance, t_snapshot, snr)

    # Matriz de covariância R = XX^H / t_snapshot
    R = np.dot(received_signal, np.conj(received_signal.T)) / t_snapshot

    # Inverter a matriz de covariância
    R_inv = np.linalg.inv(R)

    # Vetor de steering para o ângulo de interesse
    theta_target_rad = (theta_target)*(np.pi/180) # Ângulo de interesse em radianos

    
    d = 0.5  # Espaçamento normalizado entre antenas
    k = 2 * np.pi  # Número de onda normalizado
    a_target = np.exp(-1j * k * d * np.arange(m_antennas) * np.sin(theta_target_rad))
    a_target = a_target.reshape(-1, 1)

    # Cálculo dos pesos MVDR
    w_mvdr = (R_inv @ a_target) / (a_target.conj().T @ R_inv @ a_target)

    # Varrendo ângulos para construir o padrão de radiação
    theta_scan = np.linspace(-90, 90, 360)  # Ângulos de -90° a 90°
    B = np.zeros_like(theta_scan, dtype=complex)

    for i, theta in enumerate(theta_scan):
        theta_rad = np.deg2rad(theta)
        a_scan = np.exp(-1j * k * d * np.arange(m_antennas) * np.sin(theta_rad)).reshape(-1, 1)
        B[i] = w_mvdr.conj().T @ a_scan  # Beampattern

    return lin2db(np.abs(B)), theta_scan
