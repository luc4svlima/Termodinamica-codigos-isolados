#Ideia

#Dados:
#temperatura T
#composição líquida x

#Quer:
#pressão de bolha P
#composição vapor y

#Na forma ideal, nem precisa método numérico:

#𝑃 = ∑ 𝑥𝑖 𝑃𝑖𝑠𝑎𝑡(𝑇) 

# Nome: Ponto de Bolha em Pressão (Bolha P)
# O que faz? Calcula a pressão de bolha e a composição do vapor y
# Qual caminho aborda? Usa a Lei de Raoult na forma ideal

# Recomendação:
# Misturas ideais ou como base de estudo
# Pressão baixa
# Quando quero entender a estrutura do problema

# Não recomendado:
# Sistemas fortemente não ideais
# Alta pressão
# Quando preciso considerar gama ou fi

# Como funciona?
# 1 - Calcular P_sat de cada componente na temperatura dada
# 2 - Calcular P = soma(x_i * P_sat_i)
# 3 - Calcular y_i = x_i * P_sat_i / P


import numpy as np


# Antoine: ln(P_sat) = A - B / (T + C)
def psat_antoine(T, A, B, C):
    return np.exp(A - B / (T + C))


def bolha_P_ideal(T, x, A, B, C):
    x = np.array(x, dtype=float)
    A = np.array(A, dtype=float)
    B = np.array(B, dtype=float)
    C = np.array(C, dtype=float)

    P_sat = psat_antoine(T, A, B, C)

    P = np.sum(x * P_sat)

    y = x * P_sat / P

    return {
        "P_bolha": P,
        "y": y,
        "P_sat": P_sat
    }


# EXEMPLO
T = 350.0
x = [0.4, 0.6]

A = [14.2724, 14.2043]
B = [2945.47, 2972.64]
C = [224.0, 209.0]

resultado = bolha_P_ideal(T, x, A, B, C)

print("\n--- RESULTADOS ---")
print("Ponto de bolha (P):", resultado["P_bolha"])
print("Composição do vapor (y):", resultado["y"])
print("P_sat:", resultado["P_sat"])