import argparse
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from numpy import genfromtxt

def main():
    # argument parser
    parser = argparse.ArgumentParser(description="Short sample app")
    parser.add_argument("file", help="Path to the .csv file")
    args = parser.parse_args()

    # read input file data into arrays
    plotdata = genfromtxt(args.file, delimiter=',')
    x = plotdata[:, 0]
    y = plotdata[:, 1]
    
    # calculate asymptotic
    alpha0 = 6e-9;
    omegap =9.0;
    gamma  =0.1;
    rho = gamma/(omegap*omegap);
    omegaa = 1.3;
    za = 0.01;
    Gconst = alpha0*omegaa*omegaa*rho/(za*za*za);
    y2 = np.ones(len(x))*Gconst;

    # make loglogplot from this
    matplotlib.rcParams["text.usetex"] = True;
    matplotlib.rcParams["font.size"] = 24;
    # math text
    plt.figure(figsize=(10, 3))
    plt.loglog(x,y,label='decay')
    plt.loglog(x,y2,label='decay_const')
    plt.title('Decay rate over frequency')
    plt.ylabel(r'$\Gamma(\omega)$ in eV')
    plt.xlabel(r'$\omega$ in eV')
    plt.show()


if __name__ == "__main__":
    main()
