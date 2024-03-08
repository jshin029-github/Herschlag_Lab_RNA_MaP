# To CPseq file that has already been grepped for the RNAP, FID, or PX_CONTROL seqs,
# add a flag in the second column for downstream processing

import pandas as pd
import sys, getopt

# for reference, not functional
FID = "CTTGGGTCCACAGGACACTCG"
anyRNA = "TTATGCTATAATTATT"

def Add_Preset_Flag(flag, my_df):
    r,c = my_df.shape
    flag_list = []
    for i in range(r):
        flag_list.append(flag)
    #print(len(flag_list))
    #print(flag_list[:5])
    col_flag = pd.DataFrame(flag_list)
    my_df.iloc[:, 1] = col_flag.iloc[:, 0].values
    return(my_df)




#-------------main------------------

def main(argv):


    # Assign command lne args
    try:
        opts, args = getopt.getopt(argv, ":ho:i:f:", ["ifile="])
    except getopt.GetoptError:
        print('addLibFlag.py -i <inputfile>')
        sys.exit(1)
    for opt, arg in opts:
        if opt == '-h':
            print('addLibFlag.py -i <inputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputfile = arg
        elif opt in ("-f", "--flag"):
            foo_flag = arg
        elif opt in ("-o", "--ofile"):
            outputfile = arg

    print("Calling writeSetFlagtoCPseq.py on", inputfile, "with flag:", foo_flag)

    # read in file to be processed
    df = pd.read_csv(inputfile, delimiter="\t", header=None, index_col=False)

    df2 = Add_Preset_Flag(foo_flag, df)
    print(df2.iloc[:, 1])

    if outputfile:
        pass
    else:
        input = str(inputfile)
        flagy = str(foo_flag)
        outputfile = flagy + "_flagged_" + input

    df2.to_csv(outputfile, sep='\t', index=False, header = False)


if __name__ == "__main__":
    main(sys.argv[1:])
