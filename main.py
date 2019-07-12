# -*- coding:utf-8 -*-
"""Particle Swarm Optimization for Cancer Evolution

Usage:
    pso.py (--infile <infile>) [--particles <particles>] [--iterations <iterations>] [--alpha=<alpha>] [--beta=<beta>] [--k=<k>] [--c1=<c1>] [--c2=<c2>] [--seed=<seed>] [--mutfile <mutfile>] [--multiple <runptcl>...]
    pso.py -h | --help
    pso.py -v | --version

Options:
    -h --help                               Shows this screen.
    --version                               Shows version.
    -i infile --infile infile               Matrix input file.
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
import os

import numpy as np
from docopt import docopt

import matplotlib.pyplot as plt

from Data import Data
import pso
from datetime import datetime


def main(argv):
    arguments = docopt(__doc__, version="PSO-Cancer-Evolution 1.0")

    particles = int(arguments['--particles'])
    iterations = int(arguments['--iterations'])
    alpha = float(arguments['--alpha'])
    beta = float(arguments['--beta'])
    k = int(arguments['--k'])
    c1 = float(arguments['--c1'])
    c2 = float(arguments['--c2'])
    seed = int(arguments['--seed'])
    runs = list(map(int, arguments['<runptcl>']))


    with open(arguments['--infile'], 'r') as f:
        # reading the file and feeding it to numpy
        # assuring that we at least have 2D array to work with
        matrix = np.atleast_2d(np.loadtxt(io.StringIO(f.read())))

    # number of mutations = number of columns
    mutations = matrix.shape[1]
    # number of cells = number of rows
    cells = matrix.shape[0]

    if arguments['--mutfile']:
        with open(arguments['--mutfile']) as f2:
            mutation_names = [l.strip() for l in f2.readlines()]
            if len(mutation_names) != mutations:
                raise Exception("Mutation names number in file does not match mutation number in data!", len(mutation_names), mutations)
    else:
        mutation_names = [i + 1 for i in range(mutations)]
    
    if k == mutations:
        raise Exception("Cannot have same possibile losses as mutations")

    matrix = matrix.tolist()

    base_dir = "results" + datetime.now().strftime("%Y%m%d%H%M%S")

    if runs:
        runs_data = []
        for r, ptcl in enumerate(runs):
            print ("=== Run number %d ===" % r)
            run_dir = base_dir + "/p%d_i%d" % (ptcl, iterations)
            if not os.path.exists(run_dir):
                os.makedirs(run_dir)
            data, helper = pso.init(ptcl, iterations, matrix, mutations, mutation_names, cells, alpha, beta, k, c1, c2, seed)
            data.summary(helper, run_dir)
            runs_data.append(data)

        Data.runs_summary(runs, runs_data, base_dir)

    else:
        run_dir = base_dir + "/p%d_i%d" % (particles, iterations)
        if not os.path.exists(run_dir):
            os.makedirs(run_dir)
        data, helper = pso.init(particles, iterations, matrix, mutations, mutation_names, cells, alpha, beta, k, c1, c2, seed)
        data.summary(helper, run_dir)

if __name__ == "__main__":
    main(sys.argv[1:])
