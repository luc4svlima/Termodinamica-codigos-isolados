# Nome: Ponto de Bolha em Temperatura (Bolha T)
# O que faz? Calcula a temperatura de bolha
# Qual caminho aborda? Resolve F(T)=0 usando Newton-Raphson

# Recomendação:
# Misturas ideais ou como base de estudo
# Pressão conhecida
# Quando quero encontrar a temperatura de ebulição da mistura

# Não recomendado:
# Sistemas fortemente não ideais
# Mau chute inicial
# Antoine fora da faixa de validade

# Como funciona?
# 1 - Definir F(T) = soma(x_i * P_sat_i(T)) - P
# 2 - Aplicar Newton-Raphson para encontrar T
# 3 - Calcular y_i no final


import numpy as np


def psat_antoine(T, A, B, C):
    return np.exp(A - B / (T + C))


def dpsat_dT(T, A, B, C):
    P_sat = psat_antoine(T, A, B, C)
    return P_sat * (B / (T + C)**2)


def bolha_T_ideal(P, x, A, B, C, T0, tolerancia=1e-6, max_iter=100):
    x = np.array(x, dtype=float)
    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)
    C = np.array(C, dtype=float)

    T = T0

    for i in range(max_iter):
        P_sat = psat_antoine(T, A, B, C)
        dP_sat = dpsat_dT(T, A, B, C)

        f = np.sum(x * P_sat) - P
        df = np.sum(x * dP_sat)

        if df == 0:
            raise ValueError("Derivada zero. Método falhou.")

        T_novo = T - f / df
        erro = abs(T_novo - T)

        if erro < tolerancia:
            y = x * P_sat / P
            return {
                "T_bolha": T_novo,
                "y": y,
                "iteracoes": i + 1,
                "erro": erro
            }

        T = T_novo

    P_sat = psat_antoine(T, A, B, C)
    y = x * P_sat / P

    return {
        "T_bolha": T,
        "y": y,
        "iteracoes": max_iter,
        "erro": erro
    }


# EXEMPLO
P = 760.0
x = [0.4, 0.6]

A = [14.2724, 14.2043]
B = [2945.47, 2972.64]
C = [224.0, 209.0]

T0 = 350.0

resultado = bolha_T_ideal(P, x, A, B, C, T0)

print("\n--- RESULTADOS ---")
print("Temperatura de bolha (T):", resultado["T_bolha"])
print("Composição do vapor (y):", resultado["y"])
print("Número de iterações:", resultado["iteracoes"])
print("Erro final:", resultado["erro"])
