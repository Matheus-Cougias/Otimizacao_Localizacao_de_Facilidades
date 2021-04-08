import argparse
import os
from amplpy import AMPL
from amplpy import DataFrame


def main():
    args = parse_arguments()
    dat = read_datafile(args.f)
    sol = call_ampl(dat)
    print(sol[0], sol[1])
    greedy_heur(dat)


def greedy_heur(dat):
    ni = dat['ni']
    nj = dat['nj']
    c = dat['c']
    Ni = dat['a']

    I = range(ni)
    J = range(nj)
    Nj = [[i for i in I if j in Ni[i]] for j in J]
    R = set(range(ni))
    S = set()
    K = [0.0] * nj

    stop = False
    h = 0
    _S = set()
    while (stop == False):
        if len(R) > 0:
            for j in J:
                K[j] = len(set(Nj[j]).intersection(R))

            _i, _ = min([(i, len(Ni[i])) for i in R], key=lambda k: k[1])
            _j, _ = min([(j, c[j] / max(K[j], 1e-5)) for j in Ni[_i]], key=lambda k: k[1])

            R = R.difference(Nj[_j])
            S = S.union([_j])
        else:
            stop = True
            _S = get_prime_cover(I, Ni, S)
    print(sum([c[j] for j in _S]), _S)


def get_prime_cover(I, Ni, S):
    _S = set(S)
    while (len(_S) > 0):
        j = _S.pop()
        Ss = S - set([j])
        is_cover = True
        for i in I:
            if len(set(Ni[i]).intersection(Ss)) == 0:
                is_cover = False
                break
        if (is_cover == True):
            S.discard(j)
    return S


def call_ampl(dat):
    ampl = AMPL()
    ampl.read('SetCovering.md')
    ampl.setOption('solver', 'cplex')

    ni = ampl.getParameter('ni')
    ni.set(dat['ni'])

    nj = ampl.getParameter('nj')
    nj.set(dat['nj'])

    c = ampl.getParameter('c')
    c.setValues(dat['c'])
    for i in range(dat['ni']):
        ampl.set['a'][i] = dat['a'][i]

    # ampl.solve()
    ampl.getOutput('solve;')

    of = ampl.getObjective('of').value()
    x = ampl.getVariable('x')
    _x = [int(j) for j in range(dat['nj']) if x[j].value() > 0.9]

    return (of, _x)


def read_datafile(name):
    if (os.path.exists(name) == False):
        raise Exception("error reading file {:s}".format(name))
    with open(name, "r") as fl:
        lns = (ln.strip() for ln in fl)
        lns = (ln for ln in lns if ln)
        ni, nj = list(map(int, next(lns).split()))
        ic = 0
        c = []
        a = []
        while (ic < nj):
            pn = list(map(int, next(lns).split()))
            c.extend(pn)
            ic += len(pn)
        ir = 0
        while (ir < ni):
            nc = list(map(int, next(lns).split()))[0]
            ic = 0
            row = []
            while (ic < nc):
                pn = list(map(int, next(lns).split()))
                row.extend(pn)
                ic += len(pn)
            a.append(row)
            ir += 1
        for i in range(ni):
            for j in range(len(a[i])):
                a[i][j] -= 1
        dat = {}
        dat.update({'ni': ni})
        dat.update({'nj': nj})
        dat.update({'c': c})
        dat.update({'a': a})
        return dat


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('f', type=str, help='data file')
    # parser.add_argument('-n', type=int, help = 'valor adicional ')
    return parser.parse_args()


if __name__ == "__main__":
    main()