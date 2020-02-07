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

    # make loglogplot from this
    matplotlib.rcParams["text.usetex"] = True;
    matplotlib.rcParams["font.size"] = 24;
    plt.figure(figsize=(10, 3))
    plt.loglog(x,-y)
    plt.show()


if __name__ == "__main__":
    main()
