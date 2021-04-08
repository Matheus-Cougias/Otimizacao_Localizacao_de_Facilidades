"""################################# RESULTADOS ###########################"""
"""Ao realizar o teste com uma quantidade de 324 clientes, obteve-se os seguintes resultados:
    -> Instalando 5 facilidades com a distância máxima de 800: 
        AMPL = 12152
        Gulosa = 12116 (GAP = 0.2962%)
        LK = 12139 (GAP = 0.1069%) """


import argparse
import os
from math import pow
from math import sqrt
from amplpy import AMPL


def main():
    args = parse_arguments()
    dat = read_datafile(args.file1, args.file2, args.S, args.p)
    sol = call_ampl(dat)
    print('O resultado pelo AMPL foi:')
    print('A demanda total atendida foi de {}' .format(sol[0]))
    print('As facilidades instaladas foram {}' .format(sol[1]))
    print('Os clientes atendidos foram {}' .format(sol[2]))
    print('\n')
    print('O resultado pela heurística gulosa foi:')
    sol = greed_heur(dat)
    print('\n')
    print('O resultado pela heurística de Lin Kernighan foi:')
    lin_kernighan(dat, sol)


"""############################### GULOSA ############################"""
"""O funcionamento básico da heurística gulosa se dá a partir da seleção
    da melhor combinação de facilidades, selecionando uma facilidade por
    rodada até que a quantidade de facilidades selecionadas seja igual
    ao valor dado inicialmente. Então na primeira rodada ela seleciona a 
    facilidade de melhor atendimento, e nas rodadas seguintes procura
    as demais facilidades que melhor agregariam valor ao atendimento
    total."""
def greed_heur(dat):
    ni = dat['ni']
    nj = dat['nj']
    w = dat['w']
    Ni = dat['a']
    p = dat['p']

    I = range(ni)           #Recebe o tamanho do conjunto de clientes
    J = range(nj)           #Recebe o tamanho do conjunto de facilidades
    Nj = [[i for i in I if j in Ni[i]] for j in J]      #Cria a matriz com os clientes que cada facilidade consegue atender
    R = set(range(ni))      #Cria um set com os clientes
    S = set(range(nj))      #Cria um set com as facilidades
    P = set()               #Cria um set onde serão adicionadas as facilidades escolhidas
    Q = set()               #Cria um set onde serão adicionados os clientes atendidos
    K = [0.0] * nj          #Variável que calculará a demanda atendida por cada facilidade em cada rodada do while
    while len(P) < p:
        L = set(J).intersection(S)
        for j in L:
            M = set(Nj[j]).intersection(R)
            K[j] = sum(w[i] for i in M)
        _j,_ = max([(j, K[j]) for j in S], key=lambda k: k[1])
        R = R.difference(Nj[_j])            #Retira os clientes escolhidos do vetor R
        S = S.difference([_j])              #Retira a facilidade escolhida do vetor S
        P = P.union([_j])                   #Aneza a facilidade escolhida ao vetor P
        Q = Q.union(Nj[_j])                 #Anexa os clientes escolhidos ao vetor Q
    FO = sum([w[i] for i in Q])
    print('O total de demanda atendida foi de {}'.format(FO))
    print('As facilidades escolhidas foram: {}' .format(P))
    print('Os clientes atendidos foram: {}' .format(Q))

    sol = {}
    sol.update({'Q': Q})
    sol.update({'FO': FO})
    sol.update({'P': P})
    return sol

"""############################### LIN KERNIGHAN ############################"""
"""Com uma solução inicial dada pela heurística gulosa, a heurística de Lin
    Kernighan atua através do teste com outras possíveis combinações baseadas
    naquela solução inicial. Assim, ela atua trocando, a cada rodada, cada uma 
    das facilidades selecionadas inicialmente por uma nova facilidade que não 
    foi selecionada, e vai comparando qual o valor do atendimento total. Caso
    apareça uma combinação que seja melhor que a atualmente utilizada, a 
    heurística realiza a atualização dos valores, retornando ao final a melhor
    solução testada."""
def lin_kernighan (dat, sol):
    ni = dat['ni']
    nj = dat['nj']
    p = dat['p']
    w = dat['w']

    I = set(range(ni))
    J = set(range(nj))

    Ni = dat['a']
    Nj = [[i for i in I if j in Ni[i]] for j in J]

    s = dict(sol)
    _J = J.difference(s['P'])

    while len(_J) > 0:
        P = set(s['P'])
        for _k in _J:
            for _j in P:
                teste = P.difference({_j}).union({_k})
                Q = set()
                for _l in teste:
                    Q = Q.union(Nj[_l])
                FO = sum([w[i] for i in Q])
                if FO > s['FO']:
                    s['FO'] = FO
                    s['P'] = teste
                    s['Q'] = Q
            _J = _J.difference({_k})
    ss = dict(s)
    print('O total de demanda atendida foi de {}'.format(ss['FO']))
    print('As facilidades escolhidas foram: {}' .format(ss['P']))
    print('Os clientes atendidos foram: {}' .format(ss['Q']))


"""############################### FUNCAO AMPL ############################"""
"""Com os valores do .dat gerados na função abaixo, o programa já pode repassar
    as informações selecionadas para o AMPL resolver. Isso ocorre através da
    comunicação entre o Python e o AMPL, ligando este programa com o .md gerado
    de acordo com o problema. Desta maneira, ele consegue dar como retorno o 
    valor da solução ótima, os clientes atendidos e as facilidades instaladas."""
def call_ampl(dat):
    ampl = AMPL()
    ampl.read('MaxCovering.md')
    ampl.setOption('solver', 'cplex')

    ni = ampl.getParameter('ni')
    ni.set(dat['ni'])  # Envia a quantidade de clientes para o AMPL

    nj = ampl.getParameter('nj')
    nj.set(dat['nj'])  # Envia a quantidade de facilidades para o AMPL

    p = ampl.getParameter('p')
    p.set(dat['p']) # Envia o valor da quantidade máxima de facilidades para o AMPL

    w = ampl.getParameter('w')
    w.setValues(dat['w'])  # Envia os pesos de cada cliente para o AMPL

    for i in range(dat['ni']):
        ampl.set['a'][i] = dat['a'][i]  # Envia a matriz A para o AMPL

    ampl.getOutput('solve;')  # Demonstra o resultado na tela do Python

    fo = ampl.getObjective('fo').value()
    x = ampl.getVariable('x')
    _x = [int(j) for j in range(dat['nj']) if x[j].value() > 0.9]
    y = ampl.getVariable('y')
    _y = [int(i) for i in range(dat['ni']) if y[i].value() > 0.9]

    return (fo,_x, _y)


"""############################### LEITURA DOS DADOS ############################"""
"""A partir da passagem das informações iniciais, o programa busca os arquivos ditos 
    pelo usuário, onde, em caso de não encontrar, ele retorna um erro ao mesmo. Caso
    ele encontre, começa a leitura dos valores, primeiro pela quantidade de clientes
    disponíveis, e em seguida as suas localizações. Sendo que as mesmas localizações
    dos clientes são usadas para testar a implantação das facilidades, o programa
    calcula a distância euclidiana entre todos os possíveis pontos, anexando, para cada
    facilidade, os clientes que estão abaixo da distância máxima dada. Em seguida ele
    cria uma nova lista para armazenar as demandas dos clientes, de acordo com os
    valores disponíveis no segundo arquivo."""
def read_datafile(name1, name2, s, p):  # A função receberá o nome dos dois arquivos .txt além dos valores de S e de p
    if os.path.exists(name1) == False:
        raise Exception("error reading file {:s}".format(name1))
    if os.path.exists(name2) == False:
        raise Exception("error reading file {:s}".format(name2))
    with open(name1, "r") as fl:  # Começa a pesquisar dentro do primeiro arquivo
        lns = (ln.strip() for ln in fl)
        lns = (ln for ln in lns if ln)
        ni, n1, n2, n3 = list(next(lns).split())
        ni = int(ni)
        nj = ni
        PC = []
        ic = 0
        while ic < ni:
            pn = list(map(int, next(lns).split()))
            PC.append(pn)
            ic += 1
        PF = PC
        ir = 0
        a =[]
        while ir < ni:
            ik = 0
            row =[]
            while ik < nj:
                Dist = sqrt(pow((PC[ir][0] - PF[ik][0]), 2) + pow((PC[ir][1] - PF[ik][1]), 2))
                if Dist <= s:
                    pn = ik
                    row.append(pn)
                ik += 1
            a.append(row)
            ir += 1
        with open(name2, "r") as fl:  # Começa a pesquisar dentro do segundo arquivo
            lns = (ln.strip() for ln in fl)
            lns = (ln for ln in lns if ln)
            ie = 0
            w = []
            while ie < ni:
                pn = list(map(int, next(lns).split()))
                w.extend(pn)
                ie += 1
    dat = {}
    dat.update({'ni': ni})
    dat.update({'nj': nj})
    dat.update({'w': w})
    dat.update({'a': a})
    dat.update({'p': p})
    return dat


"""############################### LEITURA DA ENTRADA ############################"""
""" Nesta função, o usuário dá entrada nos nomes de:
    -> primeiro arquivo que dita as possíveis localizações a serem instaladas as facilidades
    -> segundo arquivo repassa as demandas de cada um dos pontos anteriores
    -> o terceiro valor é a distância máxima de atendimento entre uma facilidade e um cliente
    -> o quarto valor é a quantidade máxima de facilidades que serão selecionadas"""
def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('file1', type=str, help='data file')  # Faz a leitura do nome do arquivo1
    parser.add_argument('file2', type=str, help='data file')  # Faz a leitura do nome do arquivo2
    parser.add_argument('S', type=int, help='data file')  # Faz a leitura da distância máxima
    parser.add_argument('p', type=int, help='data file')  # Faz a leitura do máximo de facilidades
    return parser.parse_args()

if __name__ == "__main__":
    main()