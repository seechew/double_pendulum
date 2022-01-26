#! python3 plot_data.py -f "double_pendulum.csv" -x 1 -y 2 -k 1000

import argparse
import pandas as pd
import plotext as plt
from tqdm import tqdm

# GLOBAL DATA
x_axis = []
y_axis = []

def readcsv(filename,x_col,y_col,skip_const,chunk):

    df = pd.read_csv(filename,dtype='double',chunksize=chunk)

    for chunk in tqdm(df, desc = 'CSV chunk progress bar'):
        row = chunk[chunk.index % skip_const == 0]
        x_axis.append(row.iat[0,x_col])
        y_axis.append(row.iat[0,y_col])
    
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
        '-k','--skipdata',
        metavar="K",
        default=1,
        type=int,
        help='Skip data for plotting (default=1)')

    argparser.add_argument(
        '-c','--chunk',
        metavar="C",
        default=1,
        type=int,
        help='Reading file in chunks (default=1)')
    
    args = argparser.parse_args()
    
    readcsv(args.filename,args.x_col,args.y_col,args.skipdata,args.chunk)
    plot_me()
    
