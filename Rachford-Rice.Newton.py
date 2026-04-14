# Nome: Flash (Rachford-Rice + Newton-Raphson)
# O que faz? Encontrar a fração vapor (V) e as composições x e y
# Qual caminho aborda? Resolve F(V)=0 usando Newton-Raphson

# Recomendação:
# K conhecidos (ou já calculados)
# Função suave (Rachford-Rice é bem comportada)
# Boa estimativa inicial de V (entre 0 e 1)
# Sistema com fase vapor + líquido

# Não recomendado:
# K muito próximos de 1 (pode ficar instável)
# Sistema monofásico (todo líquido ou todo vapor)
# Mau chute inicial

# Como funciona?
# 1 - Definir z (composição global) e K
# 2 - Definir F(V)
# 3 - Aplicar Newton para encontrar V
# 4 - Calcular x e y


import numpy as np


# Função Rachford-Rice
def F(V, z, K):
    return np.sum(z * (K - 1) / (1 + V * (K - 1)))


# Derivada da função
def dF(V, z, K):
    return -np.sum(z * (K - 1)**2 / (1 + V * (K - 1))**2)


def flash_newton(z, K, V0=0.5, tolerancia=1e-6, max_iter=100):
    # z = composição global, K = coeficientes de equilíbrio

    V = V0

    for i in range(max_iter):

        f = F(V, z, K)
        df = dF(V, z, K)

        if df == 0:
            raise ValueError("Derivada zero. Método falhou.")

        V_novo = V - f / df

        erro = abs(V_novo - V)

        if erro < tolerancia:
            break

        V = V_novo

    # Cálculo das composições
    x = z / (1 + V * (K - 1))
    y = K * x

    return {
        "V": V,
        "x": x,
        "y": y,
        "iteracoes": i + 1,
        "erro": erro
    }


# ---------------------------
# EXEMPLO
# ---------------------------

# composição global
z = np.array([0.5, 0.5])

# coeficientes de equilíbrio
K = np.array([2.0, 0.5])


resultado = flash_newton(z, K)

print("\n--- RESULTADOS ---")
print("Fração vapor (V):", resultado["V"])
print("Composição líquida (x):", resultado["x"])
print("Composição vapor (y):", resultado["y"])
print("Iterações:", resultado["iteracoes"])
print("Erro:", resultado["erro"])
