import argparse
from Bio import SeqIO
import os



def get_options():

    purpose = '''This is a python script to intake a alignment file and either a text file or 
    a list of ids to remove and output an alignment without these ids

    Usage: seq_remover.py --aln <input_aln> --ids <either a .txt file or string of ids> --output <out_aln_name> '''

    parser = argparse.ArgumentParser(description=purpose,
                                     prog='seq_remover.py')

    parser.add_argument('--aln', required=True, help='Alignment file to remove ids from', type=str)
    parser.add_argument('--ids', required=True, help='Prefix of output files  (required)', type=str, nargs="+")
    parser.add_argument('--output', default="trimmed_out_aln.aln", help='Out name of new alignment', type=str)
    
    args = parser.parse_args()

    return args


def main():

    input_agrs = get_options()
    fasta_file = input_agrs.aln

    ## Get the ids to remove in the right format 
    if len(input_agrs.ids) == 1:
        first_id = input_agrs.ids[0]
        if ".txt" in first_id:
            tot_ids = []
            with open(first_id, "r") as file:
                for line in file:
                    tot_ids.append(line.strip())
        else:
            tot_ids = [first_id]
    else:
        tot_ids = input_agrs.ids

    print(tot_ids)
    if os.path.exists(input_agrs.output):
        print("Removing already present out file: ", input_agrs.output)
        os.remove(input_agrs.output)


    ## Now lets loop through the seqs and remove those with ids in tot_ids 
    with open(fasta_file) as input_handle:
        alignment = SeqIO.parse(input_handle, "fasta")
        for record in alignment:
            name = record.id
            name_grep = "^" + name + "$"
            res = any([("^" + ele + "$") == name_grep for ele in tot_ids])
            if res:
                print("Removing sequence: " + name)
            else:
                ## Append lines to the output file
                with open(input_agrs.output, "a") as aligno:
                    id_line = ">" + name + "\n"
                    seq_line = str(record.seq) + "\n"
                    aligno.write(id_line)
                    aligno.write(seq_line)


    
if __name__ == '__main__':

    main()
    print("Finished")      
    
    


