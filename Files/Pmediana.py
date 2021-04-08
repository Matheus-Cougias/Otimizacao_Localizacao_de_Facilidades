"""################################# RESULTADOS ###########################"""
"""Ao realizar o teste com uma quantidade de 100 clientes, obteve-se os seguintes resultados:
    AMPL = 5819
    Gulosa = 5891 (GAP = 1.2373%)
    LK = 5819 (GAP = 0%) """

import argparse
import os
from amplpy import AMPL
from amplpy import DataFrame

def main():
    args = parse_arguments()
    dat = read_datafile(args.file1)
    sol = call_ampl(dat)
    print('O resultado pelo AMPL foi:')
    print('O custo total mínimo foi de {}' .format(sol[0]))
    print('As facilidades instaladas foram {}' .format(sol[1]))
    print('\n')
    print('O resultado pela heurística gulosa foi:')
    sol = greedy_heur(dat)
    print('\n')
    print('O resultado pela heurística de Lin Kernighan foi:')
    lin_kernighan(dat, sol)

"""######################## LIN KERNIGHAN ######################"""
""" Utilizando a solução inicial dada pela heurística gulosa, essa
    heurística realiza a troca de cada uma das facilidades por 
    facilidades que não haviam sido selecionadas. Sendo assim, sempre
    que uma combinação de menor custo total é testada, a solução
    é atualizada e melhorada."""
def lin_kernighan(dat, sol):
    ni = dat['ni']
    nj = dat['nj']
    c = dat['c']
    p = dat['p']

    I = set(range(ni))
    J = set(range(nj))

    s = dict(sol)
    _J = J.difference(s['P'])

    H = [0.0] * ni

    while len(_J) > 0:
        P = set(s['P'])
        for _k in _J:
            for _j in P:
                FO = 0
                teste = P.difference({_j}).union({_k})
                for i in I:
                    H[i],_ = min([(j, c[i][j]) for j in teste], key=lambda k: k[1])
                    FO += c[i][H[i]]
                if FO < s['FO']:
                    s['FO'] = FO
                    s['P'] = teste
            _J = _J.difference({_k})
    ss = dict(s)
    print('O custo total foi de {}'.format(ss['FO']))
    print('As facilidades escolhidas foram: {}' .format(ss['P']))


"""############################### GULOSA ################################"""
"""Procurando uma solução inicial, a primeira facilidade escolhida é aquela
    que possui o menor somatório de distância a todos os clientes. A partir 
    dessa facilidade inicial, as demais facilidades a serem implementadas 
    são selecionadas de forma que minimizem o custo total de deslocamento
    dos clientes para a facilidade instalada mais próxima."""
def greedy_heur(dat):
    ni = dat['ni']
    nj = dat['nj']
    c = dat['c']
    p = dat['p']

    I = range(ni)
    J = range(nj)
    R = set(range(ni))      #Cria um set com as facilidades
    P = set()               #Cria um set onde serão adicionadas as facilidades escolhidas
    K = [0.0] * nj          #Variável que calculará o custo de cada facilidade inicialmente

    for j in J:
        for i in I:
            K[j] += c[i][j]
    _j, _ = min([(j, K[j]) for j in J], key=lambda k: k[1])
    R = R.difference([_j])
    P = P.union([_j])
    F = [_j] * ni           #Vetor que mostra qual facilidade atenderá qual cliente

    L = [0.0] * nj
    while len(P) < p:
        for j in R:
            G = F.copy()
            for i in I:
                if c[i][j] < c[i][G[i]]:
                    G[i] = j
                L[j] += c[i][G[i]]
        _l, _ = min([(l, L[l]) for l in R], key=lambda k: k[1])
        for i in I:
            if c[i][_l] < c[i][F[i]]:
                F[i] = _l
        R = R.difference([_l])
        P = P.union([_l])
    FO = 0
    for i in I:
        FO += c[i][F[i]]

    print('O custo total foi de {}'.format(FO))
    print('As facilidades escolhidas foram: {}' .format(P))
    sol = {}
    sol.update({'FO': FO})
    sol.update({'P': P})
    return sol


"""############################### FUNCAO AMPL ###########################"""
"""Com os valores do .dat gerados na função abaixo, o programa já pode repassar
    as informações selecionadas para o AMPL resolver. Isso ocorre através da
    comunicação entre o Python e o AMPL, ligando este programa com o .md gerado
    de acordo com o problema. Desta maneira, ele consegue dar como retorno o 
    valor da solução ótima e as facilidades instaladas"""
def call_ampl(dat):
    ampl = AMPL()
    ampl.read('Pmediana.md')
    ampl.setOption('solver', 'cplex')

    ni = ampl.getParameter('ni')
    ni.set(dat['ni'])

    nj = ampl.getParameter('nj')
    nj.set(dat['nj'])

    p = ampl.getParameter('p')
    p.set(dat['p'])

    ampl.param['c'] = {(i, j) : dat['c'][i][j] for i in range(dat['ni']) for j in range(dat['nj'])}

    ampl.getOutput('solve;')

    fo = ampl.getObjective('fo').value()
    y = ampl.getVariable('y')
    _y = [int(j) for j in range(dat['nj']) if y[j].value() > 0.9]

    return (fo,_y)

def gerar_matriz_zerada (ni, nj):
    matriz = []
    for _ in range(ni):
        matriz.append( [0] * nj )
    return matriz

def gerar_matriz (ni, nj):
    matriz = []
    for _ in range(ni):
        matriz.append( [9999] * nj )
    return matriz

def floyd_algorithim (ni, c):
    k = 0
    while k < ni:
        i = 0
        while i < ni:
            j = 0
            while j < ni:
                if c[i][k] + c[k][j] < c[i][j]:
                    c[i][j] = c[i][k] + c[k][j]
                j += 1
            i += 1
        k += 1
    return c

def distancia_zero (ni, c):
    k = 0
    while k < ni:
        c[k][k] = 0
        k += 1
    return c

"""############################# LEITURA DOS DADOS ###########################"""
""" Inicialmente, caso o arquivo exista, o programa lerá quatro três informações,
    a quantidade de clientes, o número de vértices dados e quantas facilidades
    poderão ser implementadas. Sendo assim, ele faz a leitura dos vértices e,
    para complementar as distâncias dos vértices não dados, o algoritmo de 
    Floyd é aplicado, comparando as distâncias existentes e gerando novos
    valores pras demais"""
def read_datafile(name1):
    if os.path.exists(name1) == False:
        raise Exception("error reading file {:s}".format(name1))
    with open(name1, "r") as fl:  # Começa a pesquisar dentro do primeiro arquivo
        lns = (ln.strip() for ln in fl)
        lns = (ln for ln in lns if ln)
        ni, vert, p = list(next(lns).split())
        ni = int(ni)
        vert = int(vert)
        p = int(p)
        nj = ni
        c = gerar_matriz(ni, nj)
        ic = 0
        while ic < vert:
            i, j, custo = list(map(int, next(lns).split()))
            c[i-1][j-1] = custo
            c[j-1][i-1] = custo
            ic += 1
        c = distancia_zero(ni, c)
        c = floyd_algorithim(ni, c)
    dat = {}
    dat.update({'ni': ni})
    dat.update({'nj': nj})
    dat.update({'p': p})
    dat.update({'c': c})
    return dat

"""############################### LEITURA DA ENTRADA ############################"""
""" Realiza a leitura do arquivo que passará as informações do problema."""
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('file1', type=str, help='data file')  # Faz a leitura do nome do arquivo1
    return parser.parse_args()
if __name__ == "__main__":
  main()