# -*- coding:utf-8 -*-
"""Particle Swarm Optimization for Cancer Evolution

Usage:
    pso.py (--infile <infile> -m <mutations>) [--particles <particles>] [--iterations <iterations>] [--alpha=<alpha>] [--beta=<beta>] [-k=<k>] [--mutfile <mutfile>] [--c1=<c1>] [--c2=<c2>] [--vmax=<vmax>] [--inertia=<w>]
    pso.py -h | --help
    pso.py -v | --version

Options:
    -h --help                               Shows this screen.
    --version                               Shows version.
    -i infile --infile infile               Matrix input file.
    -m mutations --mutations mutations      Number of mutations [default: 0].
    --mutfile mutfile                       Path of the mutation names. If this parameter is not used, then the mutations will be named progressively from 1 to mutations.
    -p particles --particles particles      Number of particles to use for PSO [default: 5].
    -t iterations --iterations iterations   Number of iterations [default: 3].
    --alpha=<alpha>                         False negative rate [default: 0.15].
    --beta=<beta>                           False positive rate [default: 0.001].
    -k=<k>                                  K value of Dollo(k) model used as phylogeny tree [default: 3].
    --c1=<c1>                               Learning factor 1 [default: 2].
    --c2=<c2>                               Learning factor 2 [default: 2].
    --inertia=<w>                           Inertial parameter. Affects the movement propagation given the last velocity value [default: 0.5].
    --vmax=<vmax>                           Maximum particle velocity [default: 50].
"""

from docopt import docopt
import numpy as np
import io
import sys
import pso

def main(argv):
    arguments = docopt(__doc__, version="PSO-Cancer-Evolution 1.0")

    particles = int(arguments['--particles'])
    iterations = int(arguments['--iterations'])
    mutations = int(arguments['--mutations'])
    alpha = float(arguments['--alpha'])
    beta = float(arguments['--beta'])
    k = int(arguments['-k'])
    c1 = float(arguments['--c1'])
    c2 = float(arguments['--c2'])
    inertia = float(arguments['--inertia'])
    vmax = float(arguments['--vmax'])

    with open(arguments['--infile'], 'r') as f:

        if arguments['--mutfile']:
            with open(arguments['--mutfile']) as f2:
                mutation_names = [l.strip() for l in f2.readlines()]
                if len(mutation_names) != mutations:
                    raise Exception("Mutation names number file does not match command parameter!", len(mutation_names), mutations)
        else:
            mutation_names = list(range(1, mutations))

        # reading the file and feeding it to numpy
        # assuring that we at least have 2D array to work with
        matrix = np.atleast_2d(np.loadtxt(io.StringIO(f.read())))

        if matrix.shape[1] != mutations:
            raise Exception("The number of mutations in the input file does not match command parameter!")
        
        # number of cells = number of rows
        cells = matrix.shape[0]

    pso.init(particles, iterations, matrix.tolist(), mutations, mutation_names, cells, alpha, beta, k, c1, c2, inertia, vmax)

if __name__ == "__main__":
    main(sys.argv[1:])