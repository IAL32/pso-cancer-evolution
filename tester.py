import pso
import os

particles = [5, 10, 20, 30, 40, 50, 70, 75, 100, 150, 200, 500]

for n in particles:
    os.system('mkdir results/p%d_i50 && python main.py --infile "data/scg_gawad/pat1.txt" --mutations 20 --mutfile "data/scg_gawad/pat1_mut.txt" --particles %d --iterations 50 >> results/p%d_i50/result.txt && mv results/average_particle_time.png results/average_iteration_particle_time.png results/best.gv results/best.gv.png results/p%d_i50' % (n, n, n, n))
