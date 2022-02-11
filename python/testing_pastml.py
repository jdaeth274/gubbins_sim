from pastml.acr import pastml_pipeline
import ete3
import time
import dendropy
import pandas
import argparse

def get_options():
    description = 'Find node states from cps data'
    parser = argparse.ArgumentParser(description=description,
                                     prog='testing_pastml.py')

    io_opts = parser.add_argument_group('Input/output')
    io_opts.add_argument("-g",
                         "--gubbins",
                         dest="gubbins",
                         required=True,
                         help="gubbins_loc",
                         type=str)
    io_opts.add_argument("-o",
                         "--out_prefix",
                         dest="output",
                         required=True,
                         help="location of an output directory",
                         type=str)
    io_opts.add_argument("-d", "--data",
                         dest="data",
                         help="csv of ids and states",
                         required=True,
                         type=str)
    io_opts.add_argument('-s', "--state_col",
    dest="state",default=None,
    type=str,
    help="Column name to be used as state, if none provided assume state col present")

    args = parser.parse_args()
    return (args)

def pastml_run(tree, data_loc, output):
    tic_ml = time.perf_counter()
    columns = ['state']
    pastml_pipeline(data = data_loc, columns = columns, name_column = 'state',
                     tree = tree, verbose = False,
                     out_data = output,
                    work_dir = "./", prediction_method="JOINT")
    toc_ml = time.perf_counter()

    print("ML recon took this long: %s (seconds)" % (round(toc_ml - tic_ml)))

def state_col_subset(data_loc, state_col):
    ## Load up the data and subset 
    data_df = pandas.read_csv(data_loc)
    subset_df = data_df[['id',state_col]].copy()
    subset_df.rename(columns={state_col:'state'}, inplace=True)
    state_file = "./gubbins_state_res_" + state_col + ".tsv"
    subset_df.to_csv("./gubbins_state_res_" + state_col + ".tsv",
    index=False, sep="\t",header=True)

    return state_file

def main():

    # get_input args
    input_args = get_options()
    data_loc = input_args.data
    tree_loc = input_args.gubbins
    output = input_args.output
    state_col = input_args.state

    if state_col is not None:
        data_loc = state_col_subset(data_loc, state_col)

    # recon
    pastml_run(tree_loc, data_loc, output)

if __name__ == '__main__':
    main()
    print("Done on reconstruction")
