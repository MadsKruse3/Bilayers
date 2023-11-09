import numpy as np
import os
import matplotlib.pyplot as plt

def plot_distributions(monolayer_dE_zx_distribution, monolayer_dE_zy_distribution, 
                       bilayer_dE_zx_distribution, bilayer_dE_zy_distribution):
    
    plt.figure()
    plt.scatter([xy[0] for xy in bilayer_dE_zy_distribution[:]], [xy[1] for xy in bilayer_dE_zy_distribution[:]], c='b', alpha=0.5)
    plt.axline((0, 0), slope=1, color="black", linestyle=(0, (5, 5)))
    plt.tight_layout()
    plt.savefig('anisotropy_zx_distribution.pdf')

if __name__ == "__main__":
    
    plot_distributions(monolayer_dE_zx_distribution, monolayer_dE_zy_distribution, 
                       bilayer_dE_zx_distribution, bilayer_dE_zy_distribution)

    
