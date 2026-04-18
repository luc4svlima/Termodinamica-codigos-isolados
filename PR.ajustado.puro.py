# PROPRIEDADES DO COMPOSTO PURO 
# FUNÇÃO PARA DETERMINAR PARAMETROS INICIAIS
# FUNÇÃO DE TRATAMENTO DOS PARAMETROS INICIAIS
# CALCULANDO O Z
# FUNÇÃO PARA DETERMINAR LN PHI
# FUNÇÃO FINAL PARA DETERMINAR TODOS OS DADOS
# FUNÇÃO PARA RESULTADOS
# VISUALIZAÇÃO


# === PROPRIEDADES DO COMPOSTO PURO ===

Tc = 647.1 
Pc = 220.64
omega =  0.344
T = 400
P = 10

import math

R = 0.08314  # bar.L/(mol.K)  costuma ser genérico para o sistema todo

# === FUNÇÃO PARA DETERMINAR PARAMETROS INICIAIS ===

def pr_pure_parameters(T, Tc, Pc, omega):                       
    Tr = T / Tc                                                 # temperatura reduzida
    kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega**2      # determinando k
    alpha = (1 + kappa * (1 - math.sqrt(Tr)))**2                # determinando α - usando função math.sqrt para extrair a raiz
    a = 0.45724 * (R**2 * Tc**2 / Pc) * alpha                   # determinando a para atração molecular
    b = 0.07780 * (R * Tc / Pc)                                 # determinando b para repulsão molecular 
    return {                                                    # retorno um dicionário, a partir de agora os parametros são indicados por esses
        "Tr": Tr,                                               # nomes e estão travados neles
        "kappa": kappa,
        "alpha": alpha,
        "a": a,
        "b": b,
    }

# === FUNÇÃO DE TRATAMENTO DOS PARAMETROS INICIAIS === 

def pr_AB_pure(T, P, Tc, Pc, omega):
    pure_parameters = pr_pure_parameters(T, Tc, Pc, omega)
    A = pure_parameters["a"] * P / (R**2 * T**2)
    B = pure_parameters["b"] * P / (R * T)
    return {
        "A": A,
        "B": B,
    } 


# ===== CALCULANDO O Z =====

""" 
Vamos usar o método iterativo do Newton-Rapson 
Lembre que A e B são constantes, apenas Z é variável... precisamos indicar isso na função

"""

# Estrutura básica 
def newton(f, df, x0, tolerancia=1e-6, max_iter=100):  # f = função, df = derivada, x0 = chute inicial
    
    x = x0  # Inicializa x com o chute inicial  
    
    for i in range(max_iter):  # Loop de iterações
        
        fx = f(x)  # Calcula o valor da função em x
        dfx = df(x)  # Calcula o valor da derivada em x
        
        if dfx == 0:  # Evita divisão por zero (problema crítico no método)
            raise ValueError("Derivada zero. Método falhou.")
        
        x_novo = x - fx / dfx  # Fórmula do Newton-Raphson
        
        erro = abs(x_novo - x)  # Diferença entre iterações (critério de erro)
        
        if erro < tolerancia:  # Critério de parada
            return {
                "raiz": x_novo,
                "iteracoes": i + 1,
                "erro": erro
            }
        
        x = x_novo  # Atualiza o valor para próxima iteração
    
    return {
        "raiz": x,
        "iteracoes": max_iter,
        "erro": erro
    }

# Define a função do problema para se adequar ao método
def f_PR(Z,A,B):  
    return ( Z**3 - (1- B) * Z**2 + (A - 3*B**2 - 2*B) * Z - (A*B - B**2 - B**3) )  # Z^3−(1−B)Z^2+(A−3B^2−2B)Z−(AB−B^2−B^3) = 0 


# Define a derivada da função para se adequar ao método
def df_PR(Z, A, B):
    return ( 3*Z**2 - 2*(1 - B)*Z + (A - 3*B**2 - 2*B) )  # derivada de Z^3−(1−B)Z^2+(A−3B^2−2B)Z−(AB−B^2−B^3) = 0 



# === FUNÇÃO PARA DETERMINAR LN PHI ===
"""
Cada Z, vapor ou liquido, gera um ln phi vapor e um ln phi liquido 

"""

def ln_phi_PR_puro(Z, A, B):
    termo1 = Z - 1
    termo2 = math.log(Z - B)
    termo3 = A / (2 * math.sqrt(2) * B)
    termo4 = math.log(
        (Z + (1 + math.sqrt(2)) * B) /
        (Z + (1 - math.sqrt(2)) * B)
    )

    ln_phi = termo1 - termo2 - termo3 * termo4
    return ln_phi


# === FUNÇÃO FINAL PARA DETERMINAR TODOS OS DADOS ===
"""
 Ainda não conectamos a função do método com a função da EOS, vamos usar um macete para trazer o método e ao mesmo tempo travar A e B

"""

def pr_puro(T, P, Tc, Pc, omega):

    AB = pr_AB_pure(T, P, Tc, Pc, omega)  

    Z_vapor = newton (
     lambda Z: f_PR(Z, AB["A"], AB["B"]),
     lambda Z: df_PR(Z, AB["A"], AB["B"]),
     x0=1.0                       # estamos chutando 1 para achar o Z do maior volume -> VAPOR
     )  

    Z_liquido = newton(
      lambda Z: f_PR(Z, AB["A"], AB["B"]),
      lambda Z: df_PR(Z, AB["A"], AB["B"]),
      x0=0.05                      # estamos chutando 0,05 para achar o Z do menor volume -> LIQUIDO
      )


    ln_phi_vapor = ln_phi_PR_puro(Z_vapor["raiz"], AB["A"], AB["B"])
    ln_phi_liquido = ln_phi_PR_puro(Z_liquido["raiz"], AB["A"], AB["B"])

    phi_vapor = math.exp(ln_phi_vapor)
    phi_liquido = math.exp(ln_phi_liquido)

    return {
         "Z_vapor": Z_vapor["raiz"],
         "Z_liquido": Z_liquido["raiz"],
         "phi_vapor": phi_vapor,
         "phi_liquido": phi_liquido,
         "ln_phi_vapor": ln_phi_vapor,
         "ln_phi_liquido": ln_phi_liquido,
         "A": AB["A"],         
         "B": AB["B"] 
        }

resultado = pr_puro(T, P, Tc, Pc, omega)


# === FUNÇÃO PARA RESULTADOS === 

def print_resultado(resultado):
    print("\n===== RESULTADOS PR =====")

    print(f"\nA = {resultado['A']:.6f}")
    print(f"B = {resultado['B']:.6f}")

    print(f"Z vapor     = {resultado['Z_vapor']:.6f}")
    print(f"Z líquido   = {resultado['Z_liquido']:.6f}")

    print(f"\nln(phi) vapor   = {resultado['ln_phi_vapor']:.6f}")
    print(f"ln(phi) líquido = {resultado['ln_phi_liquido']:.6f}")

    print(f"\nphi vapor   = {resultado['phi_vapor']:.6f}")
    print(f"phi líquido = {resultado['phi_liquido']:.6f}")

# === VISUALIZAÇÃO ===

print_resultado(resultado)
    
    






