# Nome: Ponto de Orvalho em Temperatura (Orvalho T)
# O que faz? Calcula a temperatura de orvalho
# Qual caminho aborda? Resolve F(T)=0 usando Newton-Raphson

# Recomendação:
# Misturas ideais ou como base de estudo
# Pressão conhecida
# Quando quero encontrar a temperatura de condensação da mistura

# Não recomendado:
# Sistemas fortemente não ideais
# Mau chute inicial
# Antoine fora da faixa de validade

# Como funciona?
# 1 - Definir F(T) = soma(y_i * P / P_sat_i(T)) - 1
# 2 - Aplicar Newton-Raphson para encontrar T
# 3 - Calcular x_i no final


import numpy as np


def psat_antoine(T, A, B, C):
    return np.exp(A - B / (T + C))


def dpsat_dT(T, A, B, C):
    P_sat = psat_antoine(T, A, B, C)
    return P_sat * (B / (T + C)**2)


def orvalho_T_ideal(P, y, A, B, C, T0, tolerancia=1e-6, max_iter=100):
    y = np.array(y, dtype=float)
    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)
    C = np.array(C, dtype=float)

    T = T0

    for i in range(max_iter):
        P_sat = psat_antoine(T, A, B, C)
        dP_sat = dpsat_dT(T, A, B, C)

        f = np.sum(y * P / P_sat) - 1.0
        df = np.sum(-y * P * dP_sat / (P_sat**2))

        if df == 0:
            raise ValueError("Derivada zero. Método falhou.")

        T_novo = T - f / df
        erro = abs(T_novo - T)

        if erro < tolerancia:
            x = y * P / P_sat
            return {
                "T_orvalho": T_novo,
                "x": x,
                "iteracoes": i + 1,
                "erro": erro
            }

        T = T_novo

    P_sat = psat_antoine(T, A, B, C)
    x = y * P / P_sat

    return {
        "T_orvalho": T,
        "x": x,
        "iteracoes": max_iter,
        "erro": erro
    }


# EXEMPLO
P = 760.0
y = [0.5, 0.5]

A = [14.2724, 14.2043]
B = [2945.47, 2972.64]
C = [224.0, 209.0]

T0 = 360.0

resultado = orvalho_T_ideal(P, y, A, B, C, T0)

print("\n--- RESULTADOS ---")
print("Temperatura de orvalho (T):", resultado["T_orvalho"])
print("Composição do líquido (x):", resultado["x"])
print("Número de iterações:", resultado["iteracoes"])
print("Erro final:", resultado["erro"])
