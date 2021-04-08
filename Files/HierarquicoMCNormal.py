"""################################# RESULTADOS ###########################"""
"""Ao realizar o teste com uma quantidade de 300 clientes, obteve-se os seguintes resultados, onde cada rodada
    significa a utilização de novos valores para demanda. Foram selecionadas 15 facilidades de baixo nível e 8
    de alto nível, onde as distâncias de atendimento são de R1=20, R2=30 e T1=25:
    -> Rodada 1: 
        AMPL = 17265
        Gulosa = 14180 (GAP = 17.8685%)
        LK = 16690 (GAP = 3.3304%) 
        
    -> Rodada 2: 
        AMPL = 16857
        Gulosa = 13789 (GAP = 16.8686%)
        LK = 16350 (GAP = 3.0076%)
        
    -> Rodada 3: 
        AMPL = 16863
        Gulosa = 14065 (GAP = 16.5925%)
        LK = 16277 (GAP = 3.4750%)"""



import argparse
import os
import numpy
import random
from amplpy import AMPL
from amplpy import DataFrame


def main():
    args = parse_arguments()
    dat = read_datafile(args.file1, args.p, args.q, args.R1, args.R2, args.T1)
    sol = call_ampl_uniforme(dat)
    print('Para o problema uniforme, o AMPL achou:')
    print('A demanda total atendida foi de {}'.format(sol[0]))
    print('As facilidades de baixo nivel instaladas foram: {}'.format(sol[2]))
    print('As facilidades de alto nivel instaladas foram: {}'.format(sol[3]))
    print('Os clientes atendidos foram: {}'.format(sol[1]))
    print('\n')
    print('A heurística gulosa achou:')
    sol = greedy_heur_uniforme(dat)
    print('\n')
    print('A heurística de Lin Kernighan achou:')
    lin_kernighan_uniforme(dat, sol)


"""############################### LIN KERNIGHAN ############################"""
""" A partir da solução inicial dada, testes são realizados trocando as facilidades
    de alto nível que foram selecionadas, no objetivo de otimizar, inicialmente,
    os clientes atendidos em ambos os níveis destas facilidades. Em seguida, com
    as de alto nível já selecionadas, um passo idêntico ao da seleção de facilidades
    de baixo nível da heurística gulosa é executado."""
def lin_kernighan_uniforme (dat, sol):
    ni = dat['ni']
    nj = dat['nj']
    a = dat['a']
    b = dat['b']
    c = dat['c']
    p = dat['p']
    q = dat['q']
    f = dat['f2']

    I = set(range(ni))
    J = set(range(nj))

    s = dict(sol)
    _J = J.difference(s['Q'])

    FuncaoObjetiva = 0
    CFatendidos = set()
    CFpossiveis = []
    for _ in range(ni):
        CFpossiveis.append([0] * nj)
    Qfinal = set()

    while len(_J) > 0:
        Q = set(s['Q'])
        for _k in _J:
            for _j in Q:
                teste = Q.difference({_k}).union({_j})
                CFat = set()
                CFp = []
                for _ in range(ni):
                    CFp.append([0] * nj)
                for i in I:
                    for j in teste:
                        if b[i][j] == 1:
                            CFat = CFat.union({i})
                        if b[i][j] == 0 and c[i][j] == 1:
                            CFp[i][j] = 1
                FO = sum([f[i] for i in CFat])
                if FO > FuncaoObjetiva:
                    FuncaoObjetiva = FO
                    CFatendidos = CFat
                    CFpossiveis = CFp
                    Qfinal = Q
            _J = _J.difference({_k})

    CopiaAtendidos = set()
    CopiaAtendidos = CopiaAtendidos.union(CFatendidos)
    P = set()
    S2 = set(range(nj))
    R = set(range(ni))

    while len(P) < p:
        Z = [0.0] * nj
        for j in S2:
            for k in Qfinal:
                Z[j] = sum(f[i]*a[i][j]*CFpossiveis[i][k] for i in R)     # Calcula a demanda atendida dentro da interseção de a e c
        _z, _ = max([(j, Z[j]) for j in S2], key=lambda k: k[1]) # Seleciona a facilidade que mais agregará clientes finais
        S2 = S2.difference([_z])
        P = P.union([_z])
        for i in R:
            for k in Qfinal:
                if a[i][_z] == 1 and CFpossiveis[i][k] == 1:
                    CFatendidos = CFatendidos.union([i])  # Anexa os clientes finais atendidos
                    CFpossiveis[i][k] = 0 # Zera, pois o cliente já será atendido
    FO = sum([f[i] for i in CFatendidos])
    FuncaoObjetiva = FO

    _J = J.difference(P)
    Pfinal = set(P)

    while len(_J) > 0:
        Pteste = set(Pfinal)
        for _k in _J:
            for _j in Pteste:
                teste = Pteste.difference({_k}).union({_j})
                CFat = set()
                for i in I:
                    for j in teste:
                        if a[i][j] == 1:
                            CFat = CFat.union({i})
                CFteste = set()
                CFteste = CFat.union(CopiaAtendidos)
                FO = sum([f[i] for i in CFteste])
                if FO > FuncaoObjetiva:
                    FuncaoObjetiva = FO
                    CFatendidos = CFteste
                    P = numpy.copy(Pteste)
            _J = _J.difference({_k})

    s['FO'] = FuncaoObjetiva
    s['Q'] = Qfinal
    s['P'] = P
    s['C'] = CFatendidos
    print('Demanda total atendida de {}'.format(s['FO']))
    print('As facilidades de baixo nível instaladas foram: {}'.format(s['P']))
    print('As facilidades de alto nível instaladas foram: {}'.format(s['Q']))
    print('Os clientes atendidos foram: {}'.format(s['C']))
    return sol


"""############################### GULOSA ############################"""
""" O primeiro passo da heurística gulosa é a seleção das melhores facilidades
    de alto nível, já que somente elas conseguem atender a este nível dos clientes.
    Cada facilidade é selecionada de modo que aumente a demanda total atendida
    pelo conjunto. Então, com estas facilidades selecionadas, e considerando que
    elas também conseguem atender no baixo nível, um conjunto de facilidades de 
    baixo nível é previamente criado com as facilidades de alto nível. Para 
    selecionar as facilidade de baixo nível, procura-se as facilidades que 
    consigam aumentar o raio de atendimento para ambos os níveis."""
def greedy_heur_uniforme(dat):
    ni = dat['ni']
    nj = dat['nj']
    a = dat['a']
    b = dat['b']
    c = dat['c']
    p = dat['p']
    q = dat['q']
    f = dat['f2']

    I = range(ni)
    J = range(nj)
    R = set(range(ni))      #Cria um set com os clientes
    S = set(range(nj))      #Cria um set com as facilidades de segundo nivel
    P = set()               #Cria um set onde serão adicionadas as facilidades de baixo nível
    Q = set()               #Cria um set onde serão adicionadas as facilidades de alto nível
    C = set()               #Cria um set onde serão adicionados os clientes atendidos
    CF = set()
    K = [0.0] * nj          #Variável que calculará a demanda atendida por cada facilidade em cada rodada do while

    while len(Q) < q:
        for j in S:
            K[j] = sum(f[i]*b[i][j] for i in R)
        _j, _ = max([(j, K[j]) for j in S], key=lambda k: k[1])
        S = S.difference([_j])
        Q = Q.union([_j])
        for i in R:
            if c[i][_j] == 1:
                R = R.difference([i])
                C = C.union([i])
            if b[i][_j] == 1:
                CF = CF.union([i])

    while len(P) < p:
        Z = [0.0] * nj
        for j in C:
            Z[j] = sum(f[i]*a[i][j] for i in R)
        _z, _ = max([(j, Z[j]) for j in C], key=lambda k: k[1])
        C = C.difference([_z])
        P = P.union([_z])
        for i in R:
            if a[i][_z] == 1:
                CF = CF.union([i])

    FO = sum([f[i] for i in CF])
    sol = {}
    sol.update({'FO': FO})
    sol.update({'P': P})
    sol.update({'Q': Q})
    sol.update({'C': CF})
    print('Demanda total atendida de {}'.format(FO))
    print('As facilidades de baixo nível instaladas foram: {}'.format(P))
    print('As facilidades de alto nível instaladas foram: {}'.format(Q))
    print('Os clientes atendidos foram: {}'.format(CF))
    return sol


"""############################### FUNCAO AMPL ############################"""
"""Com os valores do .dat gerados na função abaixo, o programa já pode repassar
    as informações selecionadas para o AMPL resolver. Isso ocorre através da
    comunicação entre o Python e o AMPL, ligando este programa com o .md gerado
    de acordo com o problema. Desta maneira, ele consegue dar como retorno o 
    valor da solução ótima, os clientes atendidos, além das facilidades  de
    baixo e alto nível instaladas."""
def call_ampl_uniforme(dat):
    ampl = AMPL()
    ampl.read('HierarquicoMC.md')
    ampl.setOption('solver', 'cplex')

    ni = ampl.getParameter('ni')
    ni.set(dat['ni'])

    nj = ampl.getParameter('nj')
    nj.set(dat['nj'])

    p = ampl.getParameter('p')
    p.set(dat['p'])

    q = ampl.getParameter('q')
    q.set(dat['q'])

    f = ampl.getParameter('f')
    f.setValues(dat['f2'])

    ampl.param['a'] = {(i, j) : dat['a'][i][j] for i in range(dat['ni']) for j in range(dat['nj'])}

    ampl.param['b'] = {(i, j) : dat['b'][i][j] for i in range(dat['ni']) for j in range(dat['nj'])}

    ampl.param['c'] = {(i, j) : dat['c'][i][j] for i in range(dat['ni']) for j in range(dat['nj'])}

    ampl.getOutput('solve;')

    fo = ampl.getObjective('fo').value()
    x = ampl.getVariable('x')
    _x = [int(i) for i in range(dat['ni']) if x[i].value() > 0.9]
    y = ampl.getVariable('y')
    _y = [int(j) for j in range(dat['nj']) if y[j].value() > 0.9]
    z = ampl.getVariable('z')
    _z = [int(j) for j in range(dat['nj']) if z[j].value() > 0.9]
    return (fo, _x, _y, _z)

def gerar_matriz (ni, nj):
    matriz = []
    for _ in range(ni):
        matriz.append( [9999] * nj )
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
    cliente são geradas a partir da distribuição estatística Normal(80,15)"""
def read_datafile(file1, p, q, R1, R2, T1):
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
        ic = 0
        while ic < vert:
            i, j, custo = list(map(int, next(lns).split()))
            d[i-1][j-1] = custo
            d[j-1][i-1] = custo
            ic += 1
        d = distancia_zero(ni, d)
        d = floyd_algorithim(ni, d)
        a = numpy.copy(d)
        b = numpy.copy(d)
        c = numpy.copy(d)
        i = 0
        while i < ni:
            j = 0
            while j < nj:
                if a[i][j] <= R1:
                    a[i][j] = 1
                else:
                    a[i][j] = 0

                if b[i][j] <= T1:
                    b[i][j] = 1
                else:
                    b[i][j] = 0

                if c[i][j] <= R2:
                    c[i][j] = 1
                else:
                    c[i][j] = 0
                j += 1
            i += 1
    f2 = []
    contador = 0
    while contador < ni:
        k = random.normalvariate(80, 15)
        f2.append(k)
        contador += 1
    dat = {}
    dat.update({'ni': ni})
    dat.update({'nj': nj})
    dat.update({'f2': f2})
    dat.update({'a': a})
    dat.update({'b': b})
    dat.update({'c': c})
    dat.update({'p': p})
    dat.update({'q': q})

    return dat


"""############################### LEITURA DA ENTRADA ############################"""
""" Faz a leitura das informações iniciais do problema, sendo elas:
    -> Arquivo que passará informações do número de clientes e os custos de atendimento
    -> Quantidade de facilidades de baixo nível a serem instaladas
    -> Quantidade de facilidades de alto nível a serem instaladas
    -> Distâncias para que o atendimento seja possível"""
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('file1', type=str, help='data file')  # Faz a leitura do nome do arquivo1
    parser.add_argument('p', type=int, help='data file')  # Faz a leitura da distância máxima
    parser.add_argument('q', type=int, help='data file')  # Faz a leitura do máximo de facilidades
    parser.add_argument('R1', type=int, help='data file')  # Faz a leitura da distância máxima
    parser.add_argument('R2', type=int, help='data file')  # Faz a leitura do máximo de facilidades
    parser.add_argument('T1', type=int, help='data file')  # Faz a leitura da distância máxima
    return parser.parse_args()

if __name__ == "__main__":
    main()