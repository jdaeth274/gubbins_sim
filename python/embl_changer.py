import os
import re
import pandas
import numpy
import time
import sys
import tqdm
import multiprocessing as mp
import argparse
import math
import psutil
import datetime

def get_options():

    purpose = '''This is a python script to intake a the embl branch recon file and output
    a 5 column table of: ['start_node','end_node','start_base','end_base','base_number'] .
    Usage: changing_embl_branch_recon_to_table.py <input_embl> <num_threads> <output_csv> '''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='Changing_embl_branch_recon_to_table.py')

    parser.add_argument('--embl', required=True, help='EMBL file to switch over', type=str)
    parser.add_argument('--output', required=True, help='Prefix of output files  (required)', type=str)
    parser.add_argument('--threads', default=1, help='Number of threads to use for the switch (must be a multiple of 6)', type=int)
    parser.add_argument('--print-file',default="./printer_output",
                        help='File to print memory debugging too',type=str)

    args = parser.parse_args()

    return args

def line_num(input_file):
    branches_and_bases = 0
    num_lines = 0
    with open(input_file, "r") as input_handle:
        for line in input_handle:
            num_lines += 1
            if re.search("variation",line):
                branches_and_bases += 1
                # base_number = re.split("variation",line)
                # base_number_only = base_number[1].rstrip()
                # base_number_only = re.sub(" ","",base_number_only)
                # base_number_only = int(base_number_only)

    return branches_and_bases, num_lines


def table_getter(input_file, num_lines, branches_and_bases, row_start, row_end, thread):
    start_array = numpy.empty((branches_and_bases, 5))
    start_array[:] = numpy.NaN
    new_table = pandas.DataFrame(data=start_array,
                                 columns=['start_node', 'end_node', 'start_base', 'end_base',
                                            'base_number'])
    # #time.sleep(thread * 5)
    # print("")
    # print(os.getpid())
    # print("On this thread: " + str(thread))
    # print(row_start)
    # print(row_end)
    # print(branches_and_bases)
    # print("##################################################")
    counter = 0
    with open(input_file, "r") as embl_file:
        for i, line in enumerate(tqdm.tqdm(embl_file, total = num_lines)) :
        #for i, line in enumerate(embl_file):
            #print("")
            #start_time = time.time()
            if i < row_start:
                continue
            elif i > row_end:
                return new_table

            # if i < row_start + 10:
            #     print("###")
            #     print(i)
            #     print(counter)
            #     flaggy = True
            # elif i % 100000 == 0:
            #     print("###")
            #     print(i)
            #     print(counter)
            # elif i > row_end - 20:
            #     print("###")
            #     print(i)
            #     print(counter)


            if re.search("variation", line):
                base_number = re.split("variation", line)
                base_number_only = base_number[1].rstrip()
                base_number_only = re.sub(" ", "", base_number_only)
                base_number_only = int(base_number_only)
                new_table.iloc[counter, 4] = base_number_only

            elif re.search("node", line):
                nodes_only = re.split("=\"",line)
                indiv_nodes = re.split("->",nodes_only[1])
                indiv_nodes[1] = indiv_nodes[1].rstrip()
                indiv_nodes[1] = indiv_nodes[1][:-1]
                indiv_nodes = pandas.Series(indiv_nodes, index=['start_node','end_node'])
                new_table.iloc[counter, [0, 1]] = indiv_nodes
                # if flaggy:
                #     print(indiv_nodes)

            elif re.search("parent_base", line):
                parents_only = re.split("=\"", line)
                indiv_base = parents_only[1]
                indiv_base = indiv_base.rstrip()
                indiv_base = indiv_base[:-1]
                #indiv_base = pandas.Series(data = str(indiv_base), index=['start_base'])
                new_table.iloc[counter, 2] = indiv_base

            elif re.search("replace", line):
                new_only = re.split("=\"", line)
                replace_base = new_only[1]
                replace_base = replace_base.rstrip()
                replace_base = replace_base[:-1]
                #replace_base = pandas.Series(data = replace_base, index=['end_base'])
                new_table.iloc[counter, 3] = replace_base
                counter += 1
                # end_time = time.time()
                # duration = end_time - start_time
                # time_left = duration * (branches_and_bases - counter)
                # if counter > 0:
                #     if counter % 100 == 0:
                #         print((counter / branches_and_bases) * 100, "%")
                #         print("This long left:", time_left)
                #         new_table.to_csv(files_for_input[2],
                #                          index=False)
            flaggy = False
    return new_table


def mp_list(num_lines, threads):
    if threads > num_lines:
        print("Number of threads greater than lines in embl file, code not set up for this yet")
        sys.exit(1)
    per_process_file_lines = math.ceil(num_lines / threads)

    # get the list of list of the finish and end lines and the start counters then iterate
    # through this list of lists
    start_lines = [0]  # start file of file in loop
    end_lines = [per_process_file_lines - 1]  # end file number to go up to
    start_counter = [math.ceil((end_lines[0] + 1) / 6)]  # the number of rows of the outfile to produce

    thread_list = [1]
    tot_list = [[start_lines[0], end_lines[0], start_counter[0], thread_list[0]]]
    for k in range(threads):
        thread_num = k + 2
        thread_list.append(thread_num)
        if thread_num > threads:
            break
        start_lines.append((per_process_file_lines * (thread_num - 1)))
        if thread_num != threads:
            end_lines.append(end_lines[len(end_lines) - 1] + per_process_file_lines)
            #current_start = (math.floor(start_lines[len(start_lines) - 1] / 6))
            sc_current = (math.ceil((end_lines[0] + 1) / 6) + 2)
            start_counter.append(sc_current)
            # Will need to add code in here for the special case when we somehow have more
            # threads than
            # if (end_lines[-1] + 1) % 6:
            #     start_counter.append()
        else:
            end_lines.append(num_lines)
            start_counter.append(math.ceil((end_lines[-1] - start_lines[-1]) / 6) + 2)

        tot_list.append([start_lines[-1], end_lines[-1], start_counter[-1], thread_list[-1]])


    return tot_list
    ## This gets out the list of lists in tot_lists for all the input for the table_getter_function


def df_combiner(df_res):
    ## Function to combine all the dfs from the cluster
    tot_df = None
    num_df = len(df_res)
    for i,df in enumerate(df_res):
        if i == 0 or tot_df is None:
            df_copy = df.copy()
        else:
            df_copy = tot_df.copy()
        df_dim = df_copy.shape
        last_row_nans = df_copy.iloc[(df_dim[0] - 1)].isnull().values.any()
        print(df_dim[0])
        if last_row_nans:
            missing = list(numpy.where(df_copy.iloc[(df_dim[0] -1)].isnull())[0])
            next_df = df_res[(i + 1)].copy()
            df_copy.iloc[(df_dim[0] - 1), missing] = next_df.iloc[0, missing]
            next_df = next_df.drop(0, axis = 0)
            tot_df = df_copy.append(next_df, ignore_index = True)
        else:
            if (i + 2) <= num_df:
                next_df = df_res[(i + 1)].copy()
                tot_df = df_copy.append(next_df, ignore_index=True)

    return tot_df






if __name__ == '__main__':

    # get the arguments

    input_args = get_options()

    # Get the line num and the number of table rows

    branches_and_bases, num_lines = line_num(input_args.embl)

    # Check if multiprocessing
    print("Combining rows into csv")
    print("")
    if input_args.threads == 1:
        newer_table = table_getter(input_args.embl, num_lines, branches_and_bases, 0, num_lines, 1)
    else:
        tot_list = mp_list(num_lines,input_args.threads)
        ## Set up the pool
        mp_pool = mp.Pool(processes=input_args.threads)
        rough_num = math.ceil(num_lines / input_args.threads)
        print_file = open(input_args.print_file, "a")
        print_file.write("Starting apply async runs " + str(datetime.datetime.now()) + "\n")
        print_file.write("Starting mem usage (MB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2) + "\n")
        print_file.close()
        df_list = [mp_pool.apply_async(table_getter, args=(input_args.embl, rough_num, rn, fs, fe,th)) for fs,fe,rn,th in tot_list]
        print_file = open(input_args.print_file, "a")
        print_file.write("Finishing apply async runs " + str(datetime.datetime.now()) + "\n")
        print_file.write("Starting mem usage (MB): " + str(psutil.Process(os.getpid()).memory_info().rss / 1024 ** 2) + "\n")
        print_file.write("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~" + "\n")
        print_file.close()
        df_res = [p.get() for p in df_list]
        df_res = [p.dropna(how = "all") for p in df_res]
        print("")
        print("")
        newer_table = df_combiner(df_res)

    print("")
    print("")
    print("Writing out results")
    newer_table.to_csv(input_args.output + ".csv",
                 index=False)
    print("Done")





