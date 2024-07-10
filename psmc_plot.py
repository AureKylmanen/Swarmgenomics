#!/usr/bin/python

#https://gtpb.github.io/PGDH19/pages/ms-psmc_practical
from matplotlib import pyplot as plt
import sys

#exampe: python psmc_plot.py name.png *.psmc

# Bin size used to generate the imput of PSMC (default is 100)
BIN_SIZE = 100

# Mutation rate per base per generation
MUTATION_RATE = 2.5e-8

# Number of years per generation
GENERAITON_TIME = 25

# Size of the plot
X_MIN = 1e3
X_MAX = 1e7
Y_MIN = 0
Y_MAX = 5e4

PLOT_PSMC_RESULTS = sys.argv


def psmc2fun(filename, s=BIN_SIZE, u=MUTATION_RATE):
    
    a = open(filename, 'r')
    result = a.read()
    a.close()

    # getting the time windows and the lambda values
    last_block = result.split('//\n')[-2]
    last_block = last_block.split('\n')
    time_windows = []
    estimated_lambdas = []
    for line in last_block:
        if line[:2]=='RS':
            time_windows.append(float(line.split('\t')[2]))
            estimated_lambdas.append(float(line.split('\t')[3]))


    # getting the estimations of theta for computing N0
    result = result.split('PA\t') # The 'PA' lines contain the values of the
                                  # estimated parameters
    result = result[-1].split('\n')[0]
    result = result.split(' ')
    theta = float(result[1])
    N0 = theta/(4*u)/s

    # Scalling times and sizes
    times = [GENERAITON_TIME * 2 * N0 * i for i in time_windows]
    sizes = [N0 * i for i in estimated_lambdas]
    
    return(times, sizes)

if __name__ == "__main__":
    
    fig = plt.figure()
    ax = fig.add_subplot(111)

    for chr in PLOT_PSMC_RESULTS[2:]:
        (estimated_times, estimated_sizes) = psmc2fun(chr, BIN_SIZE, MUTATION_RATE)    
        ax.step(estimated_times, estimated_sizes, where='post', linestyle='-', label = chr)
    
    ax.set_xlabel("Time in years")
    ax.set_ylabel("Effective size (x 10^4)")
    ax.ticklabel_format(axis='y', style='sci', scilimits=(-2,2))
    ax.grid(True)
    ax.set_xlim(X_MIN, X_MAX)
    ax.set_ylim(Y_MIN, Y_MAX)
    ax.set_xscale('log')
    plt.legend(loc = 'upper center', bbox_to_anchor=(1, 0.5))
    plt.savefig(sys.argv[1])
