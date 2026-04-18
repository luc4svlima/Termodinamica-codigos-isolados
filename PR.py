"""
Peng-Robinson (EOS que descreve o estado da substancia) para substancia pura 
O que ela faz? Calcula o fator de compressibilizade Z -> em conjunto com outros parametros -> determinamos ln phi 

Para cada componente, precisamos de:
temperatura critica 
pressão critica 
fator acêntrico      Ajuda a corrigir o comportamento real da substancia 
temperatura do sistema 
pressão do sistema 
b = representa mais o tamanho molecular e a repulsão  
a = representa mais a atração molecular  

A forma da EOS: P = RT/V−b​ − aα/V(V+b)+b(V−b)    Repare que eu tenho a, b e α... eles possuem suas próprias equações 
Sabemos que V = zRT/P

a = 0,45724 R^2 Tc^2 / Pc        b = 0,07780 R Tc/Pc   α = [1+κ​(1−(Tr)^1/2)]^2    Repare que surgiram novos parametros k e Tr 

Tr = T/Tc   k = 0,37464 + 1,54226w - 0,2699w^2 

O QUE É TABELADO ? a, b, w, Tc, Pc mas podemos ter variações 

"""

# ===== DADO UM COMPONENTE PURO, COMO O PR GERA SEUS PARÂMETROS ? =====

Tc = 647.1 
Pc = 220.64
omega =  0.344
T = 400
P = 10

import math

R = 0.08314  # bar.L/(mol.K)  costuma ser genérico para o sistema todo

def pr_pure_parameters(T, Tc, Pc, omega):                       # função para determinar o PR com base nesses parametros que eu vou definir antes
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

# ===== TRATANDO OS PARÂMETROS PARA ENCAIXAR NA EQUAÇÃO =====

"""
# A EOS não costuma ser resolvida direto com a e b 

# A = a(T) P/ R^2 T^2      B = b P/ RT

# EOS em função de z e com A e B : Z^3−(1−B)Z^2+(A−3B^2−2B)Z−(AB−B^2−B^3) = 0    Já sabemos que vai rolar um método numérico
# O que as raízes comunicam ? 
# A e B são constantes, mas fazem parte da equação. Apenas Z é variável
# Z= PV/RT     Z pequeno → V pequeno → fluido comprimido → líquido     Z grande → V grande → fluido expandido → vapor
"""

def pr_AB_pure(T, P, Tc, Pc, omega):
    pure_parameters = pr_pure_parameters(T, Tc, Pc, omega)
    A = pure_parameters["a"] * P / (R**2 * T**2)
    B = pure_parameters["b"] * P / (R * T)
    return {
        "A": A,
        "B": B,
    } 


# ===== CALCULANDO O Z =====

# Vamos usar o método iterativo do Newton-Rapson 
# Lembre que A e B são constantes, apenas Z é variável... precisamos indicar isso na função

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

# Define a função do problema
def f_PR(Z,A,B):  # só fiz a mudança para encaixar a equação cúbica 
    return ( Z**3 - (1- B) * Z**2 + (A - 3*B**2 - 2*B) * Z - (A*B - B**2 - B**3) )  # Z^3−(1−B)Z^2+(A−3B^2−2B)Z−(AB−B^2−B^3) = 0 


# Define a derivada da função
def df_PR(Z, A, B):
    return ( 3*Z**2 - 2*(1 - B)*Z + (A - 3*B**2 - 2*B) )  # derivada de Z^3−(1−B)Z^2+(A−3B^2−2B)Z−(AB−B^2−B^3) = 0 

# Ainda não conectamos a função do método com a função da EOS, vamos usar um macete para trazer o método e ao mesmo tempo
# travar A e B, considerando apenas Z como variável 

AB_pure = pr_AB_pure(T, P, Tc, Pc, omega)

Z_vapor = newton (
    lambda Z: f_PR(Z, AB_pure["A"], AB_pure["B"]),
    lambda Z: df_PR(Z, AB_pure["A"], AB_pure["B"]),
    x0=1.0                       # estamos chutando 1 para achar o Z do maior volume -> VAPOR
)  

Z_liquido = newton(
    lambda Z: f_PR(Z, AB_pure["A"], AB_pure["B"]),
    lambda Z: df_PR(Z, AB_pure["A"], AB_pure["B"]),
    x0=0.05                      # estamos chutando 0,05 para achar o Z do menor volume -> LIQUIDO
)

# ===== Determinando Ln phi =====
# Cada Z, vapor ou liquido, gera um ln phi vapor e um ln phi liquido 

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


ln_phi_vapor = ln_phi_PR_puro(Z_vapor["raiz"], AB_pure["A"], AB_pure["B"])
ln_phi_liquido = ln_phi_PR_puro(Z_liquido["raiz"], AB_pure["A"], AB_pure["B"])

phi_vapor = math.exp(ln_phi_vapor)
phi_liquido = math.exp(ln_phi_liquido)



#=====RESULTADOS EM FORMATO DE PRINT=====

#print("\nConstantes da EOS:")
#print(f"A = {AB_pure['A']:.5f}")
#print(f"B = {AB_pure['B']:.5f}")

#print("\nRaiz de vapor:")
#print(f"Z = {Z_vapor['raiz']:.4f}")
#print("Iterações =", Z_vapor["iteracoes"])
#print(f"Erro = {Z_vapor['erro']:.2e}")

#print("\nRaiz de liquido:")
#print(f"Z = {Z_liquido['raiz']:.5f}")
#print("Iterações =", Z_liquido["iteracoes"])
#print(f"Erro = {Z_liquido['erro']:.2e}")

#print("\nValores de Phi:")
#print(f"phi vapor = {phi_vapor:.5f}")
#print(f"phi liquido = {phi_liquido:.5f}")


