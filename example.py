import pandas as pd
import numpy as np
import os
print(os.getcwd())
from guide_creation import make_library


'''
Note: bowtie must be installed! Code is meant to be run from the library_design directory. 
'''

genes = ["YER018C", "YER018C"] #genes where deletions should be made
positions = [30, 90] #positions (bp) where deletions should be made
original_windows = ['GA', 'CA'] #nucleotide sequence that will be deleted 

#make_library returns a pandas dataframe with possible CRISPEY inserts for the desired deletions
#see guide_creation (and specifically make_library) for more details
inserts = make_library(genes, positions, original_windows, w=2, do_not_filter=False, deletion=True) #w=2 specifies 2-bp deletions
print(inserts)
inserts.to_csv('example_results.csv')

