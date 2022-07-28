import argparse
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

    # print force in units of micro meters / seconds^2
    hbar = 1.054571817e-34 # Joules / seconds
    c = 299792458.0 # meter / seconds
    eV = 1.602176634e-19 # Joules * seconds
    # the natural units are then
    time = hbar  / eV
    length = hbar * c / eV
    mass = eV / c**2

    # to compute the force in SI units, we have to multiply with a factor corresponding to the force in natural units
    prefactor_force = mass * length / time**2
    # go from meters to micro meters
    prefactor_force /= 1e-6

    # to get to the acceleration, we have to divide by the weight of the atom
    u = 1.66e-27 # kg
    mass_rb = 86.9 * u # weight of the rubidium atom

    # make loglogplot from this
    matplotlib.rcParams["text.usetex"] = True;
    matplotlib.rcParams["font.size"] = 24;
    plt.figure(figsize=(10, 3))
    plt.xlabel("$v / c$")
    plt.ylabel("$F$")
    plt.loglog(x,-y * prefactor_force / mass_rb)
    plt.show()


if __name__ == "__main__":
    main()
