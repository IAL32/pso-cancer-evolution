# -*- coding:utf-8 -*-
"""Particle Swarm Optimization for Cancer Evolution

Usage:
    pso.py (--infile <infile> -m <mutations>) [--particles <particles>] [--iterations <iterations>] [--alpha=<alpha>] [--beta=<beta>] [--k=<k>] [--c1=<c1>] [--c2=<c2>] [--seed=<seed>] [--mutfile <mutfile>]
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
    --beta=<beta>                           False positive rate [default: 0.00001].
    --c1=<c1>                               Learning factor for particle best [default: 0.25]
    --c2=<c2>                               Learning factor for swarm best [default: 0.75]
    --k=<k>                                 K value of Dollo(k) model used as phylogeny tree [default: 3].
    --seed=<seed>                           Seed used for RNG. If -1, a random one is generated. [default: -1]
"""

import io
import sys

import numpy as np
from docopt import docopt

import pso


def main(argv):
    arguments = docopt(__doc__, version="PSO-Cancer-Evolution 1.0")

    particles = int(arguments['--particles'])
    iterations = int(arguments['--iterations'])
    mutations = int(arguments['--mutations'])
    alpha = float(arguments['--alpha'])
    beta = float(arguments['--beta'])
    k = int(arguments['--k'])
    c1 = float(arguments['--c1'])
    c2 = float(arguments['--c2'])
    seed = int(arguments['--seed'])

    with open(arguments['--infile'], 'r') as f:

        if arguments['--mutfile']:
            with open(arguments['--mutfile']) as f2:
                mutation_names = [l.strip() for l in f2.readlines()]
                if len(mutation_names) != mutations:
                    raise Exception("Mutation names number file does not match command parameter!", len(mutation_names), mutations)
        else:
            mutation_names = [i + 1 for i in range(mutations)]

        # reading the file and feeding it to numpy
        # assuring that we at least have 2D array to work with
        matrix = np.atleast_2d(np.loadtxt(io.StringIO(f.read())))

        if matrix.shape[1] != mutations:
            raise Exception("The number of mutations in the input file does not match command parameter!")
        
        # number of cells = number of rows
        cells = matrix.shape[0]

    pso.init(particles, iterations, matrix.tolist(), mutations, mutation_names, cells, alpha, beta, k, c1, c2, seed)

if __name__ == "__main__":
    main(sys.argv[1:])
