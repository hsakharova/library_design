import pandas as pd
import numpy as np
import os
print(os.getcwd())
from guide_creation import make_library, db_genes

genes = ["YER018C", "YER018C"]
positions = [30, 90]
original_windows = ['GAG', 'GCC']


inserts = make_library(genes, positions, original_windows, w=3, do_not_filter=True, deletion=True) #delete one codon
print(inserts)

