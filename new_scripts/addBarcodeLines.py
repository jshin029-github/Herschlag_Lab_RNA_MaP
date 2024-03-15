# Written by Lena Cuevas, Stanford University, 11/18/2020
# Updates by John Shin 3/29/2022
# goal - add columns 7 and 8 to CPseq file. 7 is barcode sequecnce, 8 is barcode Phred score
# method - slices the first 16 nts from read 1 - in library 4 the first 16 nts read is my barcode
# run faster by filtering for RNAP initiation site sequence (or other lib consensus) on CPseq with grep command

# import statement
import argparse
import pandas as pd

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='Extract barcode sequence')
    parser.add_argument('--barcode_len','-b',
                        type=int, help='How many nts is barcode?')
    parser.add_argument('--CPseq_filename','-i',
                        type=str, help='Input CPseq file')
    parser.add_argument('--BC_CPseq_filename','-o',
                        type=str, help='Desired output CPseq file')

    args = parser.parse_args()
    barcode_len = args.barcode_len
    CPseq_filename = args.CPseq_filename
    BC_CPseq_filename = args.BC_CPseq_filename

    print("Calling addBarcodes.py on", CPseq_filename)

    #read in the active columns of the CPSeq to be edited
    CPseq = pd.read_csv(CPseq_filename,
                        delimiter="\t",
                        compression="gzip",
                        index_col = False,
                        usecols = [0, 1, 2, 3, 4, 5],
                        names = ["location", "exp", "R1", "QR1", "R2", "QR2"]
                        )
    # all the final column names = ["location", "exp", "R1", "QR1", "R2", "QR2", "BC", "QBC", "BS1"], )

    # CPseq.dropna(inplace = True)
    print('Input CPseq\n',CPseq.head())

    # converting to string data type
    CPseq["R1"] = CPseq["R1"].astype(str)
    CPseq["QR1"] = CPseq["QR1"].astype(str)

    # slicing away first element
    CPseq["BC"] = CPseq["R1"].str.slice(stop = (barcode_len))
    CPseq["QBC"] = CPseq["QR1"].str.slice(stop = (barcode_len))

    #checking success
    print('Output CPseq\n',CPseq.head())
    print('Size of Output CPseq',CPseq.shape)


    #make a list of dataframes and append them to the final dataframe
    BC_CPseq = CPseq[["location","exp","R1","QR1","R2","QR2","BC","QBC"]]

    #saving output file
    BC_CPseq.to_csv(BC_CPseq_filename, sep="\t", index = False, header = False)

    print('Done')
