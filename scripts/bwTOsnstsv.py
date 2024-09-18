import pyBigWig
import numpy as np
import pandas as pd
import argparse
import csv
import os

if __name__ == "__main__":
	#add a description
	parser = argparse.ArgumentParser(description="Takes Input BigWig file and outputs a Histogram Distribution of values in Chr22")

	#add the arguments
	parser.add_argument('--input_bw', dest = 'input_bw', help="input bigwig file")
	parser.add_argument('--out_name', dest = 'out_name', help="output name for files", nargs='+')

	#this allows you to access the arguments via the object args
	args = parser.parse_args()

# Open BigWig File in non-overlapping 32 bp windows

bw = pyBigWig.open(args.input_bw)
start = 0
stop = start + 32
step = 32
arr = []

# Open Values in Chr22
for i in range(50818048): #50818048

    try:
        val = bw.values("chr22", start, stop)
        arr.append(val)
    
    except:
        print("error in last")
        

    start = start + step
    stop = stop + step
    perc = (i/50818048) * 100
    print("start: ", start, "stop: ", stop, "Completed: ",perc, "% done")

arr = np.nan_to_num(arr)

# Write Array as tsv
os.makedirs(os.path.dirname(str(args.out_name[0])), exist_ok=True)

with open(str(args.out_name[0]), 'w') as f_output:
    tsv_output = csv.writer(f_output, delimiter='\t')
    tsv_output.writerows(arr)


df = pd.read_csv(str(args.out_name[0]), sep = '\t', header = None)

# Drop rows with 0s
dfdf = df.loc[~(df==0).all(axis=1)]
dfdf.reset_index(drop=True, inplace=True)

# Explode the data into columns
all_data = []
for i in range(dfdf.shape[0]):
    val = np.array(dfdf.iloc[i])
    bin_num = i
    
    baby_df = pd.DataFrame(columns=['val','bin_num'])
    baby_df['val'] = val
    baby_df['bin_num'] = bin_num
    
    all_data.append(baby_df)

DF = pd.concat(all_data, axis=0)
#DFDF = DF.loc[~(DF==0).all(axis=1)]
DFDF['log_val'] = np.log(DFDF['val'] + 1)

# Export SNS Format df for boxplot
DFDF.to_csv(str(args.out_name[1]), index =None, sep = '\t')

# Generate df intersected with bed file of GS binding sites