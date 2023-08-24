from guide_creation import *

#sorf_columns = ["chr","start","end","name","score","strand","thickstart","thickend","rgb","exons","sizes","starts"]
#sorfs = pd.read_csv('brar_data/sorfs_1.bed.txt', names=sorf_columns, delimiter='\t')
#sorfs.start = sorfs.start.astype(int)
#sorfs.end = sorfs.end.astype(int)
#sorfs.sizes = sorfs.sizes.astype(int)
#sorfs.exons = sorfs.exons.astype(int)
#print(sorfs)
#assert(np.all(sorfs.sizes == sorfs.end-sorfs.start))
#assert(np.all(sorfs.exons == 1))

#sorfs['seq'] = [get_seq(sorfs.chr[i], sorfs.start[i], sorfs.end[i]) if sorfs.strand[i]=='+'
#                else rev_complement(get_seq(sorfs.chr[i], sorfs.start[i], sorfs.end[i]))
#                for i in sorfs.index]

def validate_seq(seq):
    stop_codons = ['TAA','TAG','TGA']
    return seq[0:3]=='ATG' and seq[-3:] in stop_codons and len(seq)%3==0

def get_start_positions(seq):
    codon_seq = np.array([seq[i*3:i*3+3] for i in range(len(seq)//3)])
    start_positions = list(np.where(codon_seq == "ATG")[0])
    return start_positions #in codons

def db_genes_from_bed(db_bed):
    #ASSUMES ONLY 1 exon
    db_bed = db_bed.rename(columns={'chr':'chrom', 'start':'cdsStart', 'end':'cdsEnd'})
    if 'exons' in db_bed.columns:
        assert(np.all(db_bed.exons == 1)) #check that it's only 1
    db_bed['exonCount'] = 1
    db_bed = db_bed.set_index('name')
    assert(np.all(np.in1d(db_bed.strand, ['-','+'])))
    db_bed.strand = [1 if x=='+' else -1 for x in db_bed.strand]
    assert(np.all(db_bed.starts==0)) #can uses sizes/starts if not, but for now limit to no exons
    db_bed['exonStarts'] = [[x] for x in db_bed.cdsStart]
    db_bed['exonEnds'] = [[x] for x in db_bed.cdsEnd]
    return db_bed


#sorfs['valid_seq'] = [validate_seq(seq) for seq in sorfs.seq]
#print(np.sum(sorfs.valid_seq), "out of", len(sorfs), "are valid")
db_genes['start_positions'] = [get_start_positions(seq.replace('-','')) for seq in db_genes.sacCer3]
db_genes['num_starts'] = [len(sp) for sp in db_genes.start_positions]
print(collections.Counter(db_genes.num_starts))

#filtered for sequences that had a high number of 'N' -- 55
db_genes['num_n'] = np.array([seq.replace('-','').count('N') for seq in db_genes.sacCer3])
print('Removing db_genes with "N"')
print(collections.Counter(db_genes.num_n))
print(db_genes[db_genes.num_n != 0])
db_genes = db_genes[db_genes.num_n == 0]
#db_genes = db_genes.set_index(db_genes.name)

#ok so just -- make exons optional
#AND fix it so that 

#db_genes.to_csv('brar_data/db_genes.csv', sep='\t')



sorf_names = []
positions = []
for i in range(3): #up to 3 start codons changed
    s = db_genes[db_genes.num_starts > i]
    sorf_names += list(s.index)
    positions += [x[i]*3+1 for x in s.start_positions]#for deletion - we won't delete the A, just the TG

    #positions += [x[i]*3 for x in s.start_positions]#convert to bp
    print(f'Added {len(s.index)} targets for start number {i}')

old_windows = ['TG']*len(positions) #for deletion
#old_windows = ['ATG']*len(positions)
#new_windows = ['TAG']*len(positions)

assert(len(positions)==len(sorf_names))

print('Total targets:', len(positions))
print('Furthest start:', max(positions))
print('Setting variables')
#set_db_genes(db_genes)
#set_index_file('brar_data/sk1_original_index/sk1_o')
#set_chrom_dict(get_chrom_dict(fname='brar_data/sk1_original.fasta'))

#only the ... first 3? 
variant_inserts = make_library(sorf_names, positions, old_windows, deletion=True, w=2)

#variant_inserts = make_library(sorf_names, positions, old_windows, new_windows, deletion=False)

print('Length variants:', len(variant_inserts))
print('Variants targeting 1st position:', np.sum(variant_inserts.w_gene_pos == 1))
print('Highest position targeted (in codons)', max(variant_inserts.w_gene_pos)//3)
print('Unique genes targeted:', len(set(variant_inserts.gene)))

#variant_inserts.to_csv('brar_data/sorfs_inserts.csv')

#what db_genes is queried for?
#chrom strand exonCount exonStarts exonEnds cdsEnd
#btw and strand is -1/1


#plan:
#get deletions to work [optional - skip? leave as ''change to stop''?]
#maximize # that work
#possibility -- have the first codon deleted (frameshift)
#later codons changed to stop?


#I think for now I can just send the changed-to-stop
#unless deletions is fixed *quickly*

#then - do both sorfs and sk

#collect some info [about downstream start codons]
#try dropping start codons that are more than halfway down gene?
#(how much more if expand seed to -10? to -18?)

#ask again what their progress is like?

#hand check one? tomorrow maybe

#add neutral controls...?

#For today I'll just send the changed-to-stop versions