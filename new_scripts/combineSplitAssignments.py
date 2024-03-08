import pandas as pd

df = pd.read_csv("assign-barcodes/assignment-%d.csv"%(1))

for i in range(2,101):
    filename = "assign-barcodes/assignment-%d.csv"%(i)
    print "appending %s"%(filename)
    d = pd.read_csv(filename)
    df = df.append(d, ignore_index=True)

df.to_csv("3chips-barcode-assignments.csv", header=True, index=False)
    
