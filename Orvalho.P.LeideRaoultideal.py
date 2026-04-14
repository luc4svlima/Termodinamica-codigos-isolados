# Nome: Ponto de Orvalho em Pressão (Orvalho P)
# O que faz? Calcula a pressão de orvalho e a composição do líquido x
# Qual caminho aborda? Usa a Lei de Raoult na forma ideal

# Recomendação:
# Misturas ideais ou como base de estudo
# Pressão baixa
# Quando quero entender a condensação inicial

# Não recomendado:
# Sistemas fortemente não ideais
# Alta pressão
# Quando preciso considerar gama ou fi

# Como funciona?
# 1 - Calcular P_sat de cada componente
# 2 - Calcular 1/P = soma(y_i / P_sat_i)
# 3 - Calcular x_i = y_i * P / P_sat_i


import numpy as np


def psat_antoine(T, A, B, C):
    return np.exp(A - B / (T + C))


def orvalho_P_ideal(T, y, A, B, C):
    y = np.array(y, dtype=float)
    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)
    C = np.array(C, dtype=float)

    P_sat = psat_antoine(T, A, B, C)

    P = 1.0 / np.sum(y / P_sat)

    x = y * P / P_sat

    return {
        "P_orvalho": P,
        "x": x,
        "P_sat": P_sat
    }


# EXEMPLO
T = 350.0
y = [0.5, 0.5]

A = [14.2724, 14.2043]
B = [2945.47, 2972.64]
C = [224.0, 209.0]

resultado = orvalho_P_ideal(T, y, A, B, C)

print("\n--- RESULTADOS ---")
print("Ponto de orvalho (P):", resultado["P_orvalho"])
print("Composição do líquido (x):", resultado["x"])
print("P_sat:", resultado["P_sat"])
