#! python3 plot_data.py -f "double_pendulum.csv" -x 1 -y 2 -s

import argparse
import csv
import plotext as plt

# GLOBAL DATA
x_axis = []
y_axis = []

def readcsv(filename,skip,x_col,y_col,delim):

    with open(filename,'r') as file:
        lines = csv.reader(file,delimiter=delim)
        if skip == True: next(lines)  # SKIP HEADER

        for row in lines:
            x_axis.append(float(row[x_col]))
            y_axis.append(float(row[y_col]))

def plot_me():
    plt.colorless()
    plt.plot(x_axis,y_axis)
    plt.show()


if __name__ == "__main__":

    argparser = argparse.ArgumentParser(description='Terminal Plotting Tool')

    argparser.add_argument(
        '-f', '--filename',
        metavar='F',
        default=None,
        help='CSV filename (default=None)')

    argparser.add_argument(
        '-s', '--skip',
        action='store_true',
        help='Skip CSV file header (default=True)')

    argparser.add_argument(
        '-x', '--x_col',
        metavar="X",
        default=0,
        type=int,
        help='X-axis column number (default=0)')

    argparser.add_argument(
        '-y','--y_col',
        metavar="Y",
        default=1,
        type=int,
        help='Y-axis column number (deafult=1)')

    argparser.add_argument(
        '-d','--delim',
        metavar="D",
        default=',',
        help='Delimiter (default=",")')

    args = argparser.parse_args()
    
    readcsv(args.filename,args.skip,args.x_col,args.y_col,args.delim)
    plot_me()
    
