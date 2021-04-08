"""################################# RESULTADOS ###########################"""
"""Ao realizar o teste com uma quantidade de 300 clientes, obteve-se os seguintes resultados, onde cada rodada
    significa a utilização de novos valores para demanda. Foram selecionadas 15 facilidades de baixo nível e 8
    de alto nível, onde as distâncias de atendimento são de R1=20, R2=30 e T1=25 e a penalização=10:
    -> Rodada 1: 
        AMPL = 3348007
        Gulosa = 3977674 (GAP = 18.8072%)
        LK = 3862749 (GAP = 15.3745%) 

    -> Rodada 2: 
        AMPL = 3299265
        Gulosa = 4001446 (GAP = 21.2829%)
        LK = 3889764 (GAP = 17.8978%)

    -> Rodada 3: 
        AMPL = 3440850
        Gulosa = 4029436 (GAP = 17.0186%)
        LK = 3912171 (GAP = 13.6978%)"""




import argparse
import os
import numpy
import random
from amplpy import AMPL
from amplpy import DataFrame

def main():
    args = parse_arguments()
    dat = read_datafile(args.file1, args.p, args.q, args.R1, args.R2, args.T1, args.beta)
    print('Para o problema uniforme, o AMPL achou:')
    sol = call_ampl_uniforme(dat)
    print('O custo total foi de {}'.format(sol[0]))
    print('As facilidades de alto nivel instaladas foram: {}'.format(sol[2]))
    print('As facilidades de baixo nivel instaladas foram: {}'.format(sol[1]))
    print('\n')
    print('A heurística gulosa achou:')
    sol = greedy_heur_uniforme(dat)
    print('\n')
    print('O resultado pela heurística de Lin Kernighan foi:')
    lin_kernighan_uniforme(dat, sol)


"""############################### LIN KERNIGHAN ############################"""
def lin_kernighan_uniforme(dat, sol):
    ni = dat['ni']
    nj = dat['nj']
    p = dat['p']
    q = dat['q']
    a = dat['a']
    b = dat['b']
    c = dat['c']


    I = set(range(ni))
    J = set(range(nj))

    s = dict(sol)
    _J = J.difference(s['Q'])

    H = [0.0] * ni

    melhorQ = set(s['Q'])
    melhorFOalto = float('inf')


    while len(_J) > 0:
        Q = melhorQ
        for _k in _J:
            for _j in Q:
                FOalto = 0
                teste = Q.difference({_j}).union({_k})
                for i in I:
                    H[i],_ = min([(j, b[i][j]) for j in teste], key=lambda k: k[1])
                    FOalto += b[i][H[i]]
                if FOalto < melhorFOalto:
                    melhorFOalto = FOalto
                    melhorQ = teste
            _J = _J.difference({_k})


    custo = [float('inf')] * ni
    F = [0.0] * ni
    melhorP = set(melhorQ)
    for i in I:
        for j in melhorP:
            if c[i][j] < custo[i]:
                F[i] = j

    R = set(range(ni)).difference(melhorQ)
    while len(melhorP) < p+q:
        L = [0.0] * nj
        for j in R:
            G = F.copy()
            for i in I:
                if G[i] in melhorQ:
                    if a[i][j] < c[i][G[i]]:
                        G[i] = j
                else:
                    if a[i][j] < a[i][G[i]]:
                        G[i] = j
            for i in I:
                if G[i] in melhorQ:
                    L[j] += c[i][G[i]]
                else:
                    L[j] += a[i][G[i]]
        _l, _ = min([(l, L[l]) for l in R], key=lambda k: k[1])
        for i in I:
            if F[i] in melhorQ:
                if a[i][_l] < c[i][F[i]]:
                    F[i] = _l
            else:
                if a[i][_l] < a[i][F[i]]:
                    F[i] = _l
        R = R.difference([_l])
        melhorP = melhorP.union([_l])
    FObaixo = 0
    for i in I:
        if F[i] in melhorQ:
            FObaixo += c[i][F[i]]
        else:
            FObaixo += a[i][F[i]]

    melhorP = melhorP.difference(melhorQ)
    melhorFO = melhorFOalto + FObaixo
    _J = J.difference(melhorQ).difference(melhorP)

    custos = [float('inf')] * ni
    F = [0.0] * ni
    while len(_J) > 0:
        for _k in _J:
            for _j in melhorP:
                FO = 0
                P = melhorP.difference({_j}).union({_k}).union(melhorQ)
                for j in P:
                    for i in I:
                        soma = 0
                        if j in melhorQ:
                            if c[i][j]  < custos[i]:
                                F[i] = j
                                soma = c[i][j]
                        else:
                            if a[i][j] < custos[i]:
                                F[i] = j
                                soma = a[i][j]
                        FO += soma
                if FO < melhorFO:
                    melhorFO = FO
                    melhorP = P
            _J = _J.difference({_k})
    s['FO'] = melhorFO
    s['P'] = melhorP
    s['Q'] = melhorQ
    print('O custo total de atendimento foi de: {}' .format(s['FO']))
    print('As facilidades de alto nível escolhidas foram: {}' .format(s['Q']))
    print('As facilidade de baixo nível escolhidas foram: {}' .format(s['P']))



"""############################### GULOSA ############################"""
""" O funcionamento da heurística gulosa é simular ao da heurística gulosa
    do problema não-hierarquizado. O primeiro passo se dá na seleção da 
    facilidade de alto nível que possui menor somatorio de distâncias em 
    relação aos clientes. Assim, as demais facilidades de alto nível
    são selecionadas de modo que o custo total de deslocamento, com 
    a penalização da distância máxima de atendimento, seja minimizada.
    Sendo assim, com as facilidades de alto nível atendendo nos dois 
    níveis, um conjunto de baixo nível já é criado considerando esta
    característica do problema e as facilidades de baixo nível são
    selecionadas de modo que otimizem o atendimento para seu nível."""
def greedy_heur_uniforme(dat):
    ni = dat['ni']
    nj = dat['nj']
    p = dat['p']
    q = dat['q']
    a = dat['a']
    b = dat['b']
    c = dat['c']


    I = range(ni)
    J = range(nj)
    R = set(range(ni))
    Q = set()
    K = [0.0] * nj

    for j in J:
        for i in I:
            K[j] += b[i][j]
    _j, _ = min([(j, K[j]) for j in J], key=lambda k: k[1])
    R = R.difference([_j])
    Q = Q.union([_j])
    F = [_j] * ni

    L = [0.0] * nj
    while len(Q) < q:
        for j in R:
            G = F.copy()
            for i in I:
                if b[i][j] < b[i][G[i]]:
                    G[i] = j
                L[j] += b[i][G[i]]
        _l, _ = min([(l, L[l]) for l in R], key=lambda k: k[1])
        for i in I:
            if b[i][_l] < b[i][F[i]]:
                F[i] = _l
        R = R.difference([_l])
        Q = Q.union([_l])
    FOalto = 0
    for i in I:
        FOalto += b[i][F[i]]

    custo = [float('inf')] * ni
    F = [0.0] * ni
    P = set(Q)
    for i in I:
        for j in P:
            if c[i][j] < custo[i]:
                F[i] = j


    R = set(range(ni)).difference(Q)
    while len(P) < p+q:
        L = [0.0] * nj
        for j in R:
            G = F.copy()
            for i in I:
                if G[i] in Q:
                    if a[i][j] < b[i][G[i]]:
                        G[i] = j
                else:
                    if a[i][j] < a[i][G[i]]:
                        G[i] = j
            for i in I:
                if G[i] in Q:
                    L[j] += c[i][G[i]]
                else:
                    L[j] += a[i][G[i]]
        _l, _ = min([(l, L[l]) for l in R], key=lambda k: k[1])
        for i in I:
            if F[i] in Q:
                if a[i][_l] < c[i][F[i]]:
                    F[i] = _l
            else:
                if a[i][_l] < a[i][F[i]]:
                    F[i] = _l
        R = R.difference([_l])
        P = P.union([_l])
    FObaixo = 0
    for i in I:
        if F[i] in Q:
            FObaixo += c[i][F[i]]
        else:
            FObaixo += a[i][F[i]]
    P = P.difference(Q)
    FO = FObaixo + FOalto
    sol = {}
    sol.update({'FO': FO})
    sol.update({'FOalto' : FOalto})
    sol.update({'FObaixo' : FObaixo})
    sol.update({'Q' : Q})
    sol.update({'P': P})
    print('O custo total de atendimento foi de: {}' .format(sol['FO']))
    print('As facilidades de alto nível escolhidas foram: {}' .format(sol['Q']))
    print('As facilidade de baixo nível escolhidas foram: {}' .format(sol['P']))
    return sol



"""############################### FUNCAO AMPL ############################"""
"""Com os valores do .dat gerados na função abaixo, o programa já pode repassar
    as informações selecionadas para o AMPL resolver. Isso ocorre através da
    comunicação entre o Python e o AMPL, ligando este programa com o .md gerado
    de acordo com o problema. Desta maneira, ele consegue dar como retorno o 
    valor da solução ótima, além das facilidades  de baixo e alto nível instaladas."""
def call_ampl_uniforme(dat):
    ampl = AMPL()
    ampl.read('HierarquicoPMediana.md')
    ampl.setOption('solver', 'cplex')

    ni = ampl.getParameter('ni')
    ni.set(dat['ni'])

    nj = ampl.getParameter('nj')
    nj.set(dat['nj'])

    p = ampl.getParameter('p')
    p.set(dat['p'])

    q = ampl.getParameter('q')
    q.set(dat['q'])

    ampl.param['a'] = {(i, j) : dat['a'][i][j] for i in range(dat['ni']) for j in range(dat['nj'])}

    ampl.param['b'] = {(i, j) : dat['b'][i][j] for i in range(dat['ni']) for j in range(dat['nj'])}

    ampl.param['c'] = {(i, j) : dat['c'][i][j] for i in range(dat['ni']) for j in range(dat['nj'])}

    ampl.getOutput('solve;')

    fo = ampl.getObjective('fo').value()
    y = ampl.getVariable('y')
    _y = [int(j) for j in range(dat['nj']) if y[j].value() > 0.9]
    z = ampl.getVariable('z')
    _z = [int(j) for j in range(dat['nj']) if z[j].value() > 0.9]

    return (fo, _y, _z)



def gerar_matriz (ni, nj):
    matriz = []
    for _ in range(ni):
        matriz.append([float('inf')] * nj)
    return matriz

def floyd_algorithim (ni, d):
    k = 0
    while k < ni:
        i = 0
        while i < ni:
            j = 0
            while j < ni:
                if d[i][k] + d[k][j] < d[i][j]:
                    d[i][j] = d[i][k] + d[k][j]
                j += 1
            i += 1
        k += 1
    return d

def distancia_zero (ni, d):
    k = 0
    while k < ni:
        d[k][k] = 0
        k += 1
    return d


"""############################### LEITURA DOS DADOS ############################"""
""" Com a existência do arquivo repassado pelo usuário, os dados da quantidade de
    clientes e as distâncias entre os vértices são lidos, juntamente com as 
    distâncias disponíveis para leitura. As distâncias não dadas são geradas a 
    partir da utilização do algoritmo de Floyd, enquanto as demandas de cada
    cliente são geradas a partir da distribuição estatística Uniforme(20,30)"""
def read_datafile(file1, p, q, R1, R2, T1, beta):
    if os.path.exists(file1) == False:
        raise Exception("error reading file {:s}".format(file1))
    with open(file1, "r") as fl:
        lns = (ln.strip() for ln in fl)
        lns = (ln for ln in lns if ln)
        ni, vert, _ = list(next(lns).split())
        ni = int(ni)
        vert = int(vert)
        nj = ni
        d = gerar_matriz(ni, nj)
        a = numpy.copy(d)
        b = numpy.copy(d)
        c = numpy.copy(d)
        ic = 0
        while ic < vert:
            i, j, custo = list(map(int, next(lns).split()))
            d[i-1][j-1] = custo
            d[j-1][i-1] = custo
            ic += 1
        d = distancia_zero(ni, d)
        d = floyd_algorithim(ni, d)
        Pen1 = numpy.copy(d)
        Pen2 = numpy.copy(d)
        Pen3 =numpy.copy(d)
        while i < ni:
            j = 0
            while j < nj:
                if Pen1[i][j] <= R1:
                    Pen1[i][j] = beta
                else:
                    Pen1[i][j] = 1

                if Pen2[i][j] <= T1:
                    Pen2[i][j] = beta
                else:
                    Pen2[i][j] = 1

                if Pen3[i][j] <= R2:
                    Pen3[i][j] = beta
                else:
                    Pen3[i][j] = 1
                j += 1
            i += 1
    f = []
    contador = 0
    while contador < ni:
        j = random.uniform(20, 30)
        f.append(j)
        contador += 1
    for i in range(ni):
        for j in range(nj):
            a[i][j] = f[i] * d[i][j] * Pen1[i][j]
            b[i][j] = f[i] * d[i][j] * Pen3[i][j]
            c[i][j] = f[i] * d[i][j] * Pen2[i][j]
    dat = {}
    dat.update({'ni': ni})
    dat.update({'nj': nj})
    dat.update({'p': p})
    dat.update({'q': q})
    dat.update({'a': a})
    dat.update({'b': b})
    dat.update({'c': c})
    return dat


"""############################### LEITURA DA ENTRADA ############################"""
""" Faz a leitura das informações iniciais do problema, sendo elas:
    -> Arquivo que passará informações do número de clientes e os custos de atendimento
    -> Quantidade de facilidades de baixo nível a serem instaladas
    -> Quantidade de facilidades de alto nível a serem instaladas
    -> Distâncias para que ocorram penalizações
    -> Multiplicador de penalização"""
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('file1', type=str, help='data file')
    parser.add_argument('p', type=int, help='data file')
    parser.add_argument('q', type=int, help='data file')
    parser.add_argument('R1', type=int, help='data file')
    parser.add_argument('R2', type=int, help='data file')
    parser.add_argument('T1', type=int, help='data file')
    parser.add_argument('beta', type=int, help='data file')
    return parser.parse_args()
if __name__ == "__main__":
    main()