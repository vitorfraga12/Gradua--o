import numpy as np

# Função que transforma Linear para dB
def lin2db(x):
    return 10 * np.log10(x)

# Função que transforma dB para Linear
def db2lin(x):
    return 10 ** (x / 10)

####################################################################################################################


# Função que gera as coordenadas dos APs
def distribuir_APs(num_aps, area):
    '''Distributes Access Points (APs) evenly within a square area.
    
    Parameters:
    num_aps (int): The number of APs to distribute. Must be a perfect square.
    area (int): The length of the side of the square area in which to distribute the APs.
    
    Returns:
    np.array: An array of coordinates for the APs, or None if num_aps is not a perfect square.'''
    
    if num_aps not in [1, 4, 9, 16, 25, 36, 49, 64, 100]:
        return None

    tamanho_quadrado = area
    lado_quadrado = int(np.sqrt(num_aps))

    tamanho_celula = tamanho_quadrado // lado_quadrado

    # Criar coordenadas usando meshgrid
    x, y = np.meshgrid(np.arange(0.5 * tamanho_celula, tamanho_quadrado, tamanho_celula),
                      np.arange(0.5 * tamanho_celula, tamanho_quadrado, tamanho_celula))

    coordenadas_APs = np.column_stack((x.ravel(), y.ravel()))

    return coordenadas_APs

####################################################################################################################

# Função que gera a distância entre a UE e a AP
def dAPUE(x_coord, y_coord, ap_coord, d_reference=1):
  '''Calculate the Euclidean distance between a user equipment (UE) and an access point (AP).
    
    Parameters:
    ue_coords (tuple): A tuple (x_coord, y_coord) representing the coordinates of the UE.
    ap_coords (np.array): An array containing the coordinates of the APs.
    
    Returns:
    float: The Euclidean distance between the UE and the AP. If the euclidean distance is less than 1, return 1.
  '''
  dist = np.linalg.norm(np.array([x_coord, y_coord]) - ap_coord)
  if dist < d_reference:
    dist = d_reference

  return dist


####################################################################################################################

def find_pathgain(dist, num_channels):
    """ 
    Função que calcula o path gain entre a UE e o AP.

    Parameters:
    dist (np.array): A distância entre a UE e o AP.
    num_channels (int): O número de canais disponíveis.

    Returns:
    np.array: O path gain entre a UE e o AP.
    """


    ambiente_const = 1e-4
    pathloss_const = 4

    dist_expanded = np.repeat((dist ** pathloss_const)[:, :, np.newaxis], num_channels, axis=2)  # Agora tem forma (ues, aps, num_channels)

    # Calculando o path gain
    path_gain_result =  (ambiente_const / dist_expanded)

    return path_gain_result

####################################################################################################################

# Função que define o canal aleatório para cada UE

def ue_channel_random(distance_APUE, num_channels):
    path_gain = find_pathgain(distance_APUE, num_channels)

    # Garantindo que cada UE use um canal aleatório
    if path_gain.shape[0] > 1:
        for ue in range(path_gain.shape[0]):  # Para cada UE
            channel = np.random.randint(0, num_channels)  # Seleciona um canal aleatório
            mask = np.zeros(num_channels)
            mask[channel] = 1  # Habilita apenas o canal selecionado
            path_gain[ue, :, :] *= mask  # Aplica a máscara de canal

    return path_gain

####################################################################################################################

def calculate_sinr(banda, K_0, aps, ues, channels, p_t, area):
    '''Função que calcula a SINR para múltiplos UEs e APs em diferentes canais, 
    alocando cada UE ao AP com o maior path gain.

    Parâmetros:
    banda (float): Largura de banda total.
    K_0 (float): Constante de ruído.
    aps (int): Número de APs.
    ues (int): Número de UEs.
    channels (int): Número de canais.
    p_t (float): Potência de transmissão.
    area (float): Área de distribuição dos APs.
    
    Retorna:
    np.array: O valor do SINR para cada UE.
    '''
    # Inicializações
    x_coord, y_coord, aps_position = np.zeros(ues), np.zeros(ues), distribuir_APs(aps, area)
    dist = np.zeros((ues, aps))
    
    power_trans = np.ones(ues) * p_t
    power_noise = np.ones(ues) * (K_0 * banda / channels)

    # Coordenadas dos UEs
    for ue_index in range(ues):
        x_coord[ue_index] = np.random.randint(0, area)
        y_coord[ue_index] = np.random.randint(0, area)

    # Distância entre UEs e APs
    for ue_index in range(ues):
        for ap_index in range(aps):
            dist[ue_index][ap_index] = dAPUE(x_coord[ue_index], y_coord[ue_index], aps_position[ap_index])

    # Cálculo do path gain
    path_gain = ue_channel_random(dist, channels)

    # Alocação do AP com maior path gain para cada UE
    best_ap_index = np.argmax(np.max(path_gain, axis=2), axis=1)  # Índice do AP com maior path gain em qualquer canal

    # Cálculo do SINR
    sinr_ue = np.zeros(ues)

    for ue_index in range(ues):
        ap_index = best_ap_index[ue_index]  # Aloca o UE ao AP com maior path gain

        for channel_index in range(channels):
            # Calculando a potência recebida no AP alocado
            power_received = np.abs(path_gain[ue_index, ap_index, channel_index]) * power_trans[ue_index]

            # Calculando a interferência no mesmo canal
            interference_sum = 0
            for other_ue_index in range(ues):
                if other_ue_index != ue_index:
                    interference_sum += np.abs(path_gain[other_ue_index, ap_index, channel_index]) * power_trans[other_ue_index]

            # SINR para o canal atual
            sinr = power_received / (interference_sum + power_noise[ue_index])

            # Mantém o maior SINR entre os canais
            sinr_ue[ue_index] = max(sinr_ue[ue_index], sinr)

    return sinr_ue

####################################################################################################################

def find_capacity(sinr, banda, channels):
    '''Função que calcula a capacidade de um dado canal.
    
    Parâmetros:
    sinr (list): A relação sinal ruído mais interferência.
    banda (int): A largura de banda do canal.    
    channels (int): O número de canais.
    Retorna:
    list: A capacidade do canal.'''

    banda_canal = banda / channels  
    
    capacity = banda_canal * np.log2(1 + sinr)
    
    return capacity