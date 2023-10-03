from guide_creation import *

#TODO - list # of starts downstream from each position
#TODO - list % of gene that position is at
#TODO - filter lists, choosing which positions to edit based on - replicates, 
#TODO - add neutral edits? [but sometimes they aren't neutral...?]

set_index_file('brar_data/sk1_original_index/sk1_o')
set_chrom_dict(get_chrom_dict(fname='brar_data/sk1_original.fasta'))
#print(index_file) 

sorf_columns = ["chr","start","end","name","score","strand","thickstart","thickend","rgb","exons","sizes","starts"]
sorfs = pd.read_csv('brar_data/sorfs_1.bed.txt', names=sorf_columns, delimiter='\t')
sorfs.start = sorfs.start.astype(int)
sorfs.end = sorfs.end.astype(int)
sorfs.sizes = sorfs.sizes.astype(int)
sorfs.exons = sorfs.exons.astype(int)
#print(sorfs)
assert(np.all(sorfs.sizes == sorfs.end-sorfs.start))
assert(np.all(sorfs.exons == 1))

sorfs['seq'] = [get_seq(sorfs.chr[i], sorfs.start[i], sorfs.end[i]) if sorfs.strand[i]=='+'
                else rev_complement(get_seq(sorfs.chr[i], sorfs.start[i], sorfs.end[i]))
                for i in sorfs.index]

def validate_seq(seq):
    stop_codons = ['TAA','TAG','TGA']
    return seq[0:3]=='ATG' and seq[-3:] in stop_codons and len(seq)%3==0



def db_genes_from_bed(db_bed):
    db_bed = db_bed.rename(columns={'chr':'chrom', 'start':'cdsStart', 'end':'cdsEnd', 'exons':'exonCount'})
    db_bed = db_bed.set_index('name')
    assert(np.all(np.in1d(db_bed.strand, ['-','+'])))
    db_bed.strand = [1 if x=='+' else -1 for x in db_bed.strand]
    assert(np.all(db_bed.starts==0)) #can uses sizes/starts if not, but for now limit to no exons
    db_bed['exonStarts'] = [[x] for x in db_bed.cdsStart]
    db_bed['exonEnds'] = [[x] for x in db_bed.cdsEnd]
    return db_bed


sorfs['valid_seq'] = [validate_seq(seq) for seq in sorfs.seq]
print(np.sum(sorfs.valid_seq), "out of", len(sorfs), "are valid")
sorfs['start_positions'] = [get_start_positions(seq) for seq in sorfs.seq]
sorfs['num_starts'] = [len(sp) for sp in sorfs.start_positions]
print('number of genes with x starts:', collections.Counter(sorfs.num_starts))

#filtered for sequences that had a high number of 'N' -- 55
sorfs['num_n'] = np.array([seq.count('N') for seq in sorfs.seq])
print('Removing sorfs with "N"')
print('numbers of n:', collections.Counter(sorfs.num_n))
print(sorfs[sorfs.num_n != 0].index)
sorfs = sorfs[sorfs.num_n == 0]
sorfs = sorfs.set_index(sorfs.name)

#ok so just -- make exons optional
#AND fix it so that 

sorfs.to_csv('brar_data/sorfs.csv', sep='\t')

#now just need to add 

db_sorfs = db_genes_from_bed(sorfs)
db_sorfs.to_csv('brar_data/db_sorfs.csv')

sorf_names = []
positions = []
for i in range(11): #change to 10 ---> this is all starts *except* there was one gene that had 27 starts that I ignored
    s = sorfs[sorfs.num_starts > i]
    sorf_names += list(s.index)
    positions += [x[i]*3 for x in s.start_positions]
    #positions += [x[i]*3+1 for x in s.start_positions]#for deletion - we won't delete the A, just the TG

    #positions += [x[i]*3 for x in s.start_positions]#convert to bp
    print(f'Added {len(s.index)} targets for start number {i}')

#old_windows = ['TG']*len(positions) #for deletion
old_windows = ['ATG']*len(positions)
new_windows = ['TAA']*len(positions) #more disruptions this way?

assert(len(positions)==len(sorf_names))

print('Total targets:', len(positions))
print('Furthest start:', max(positions))
print('Setting variables')
set_db_genes(db_sorfs)

#check length 
sorf_sizes = db_sorfs.cdsEnd - db_sorfs.cdsStart
print(f"Sizes of sorfs: min{min(sorf_sizes)//3}, max {max(sorf_sizes)//3},\
      mean {np.mean(sorf_sizes)}, median {np.median(sorf_sizes)}, 10% {np.quantile(sorf_sizes, 0.1)}")
#let's see how many have guides to begin with
guide_list = []
for s in db_sorfs.index:
    #guide has to be in first half, I think
    guide_list += [(s,) + x for x in 
                   get_all_guides(db_sorfs.chrom[s], db_sorfs.cdsStart[s], end=db_sorfs.cdsEnd[s]-3)] 
guides_db = pd.DataFrame(guide_list, 
                              columns=["gene","guide_seq","chrom","cut_pos","guide_strand"])
guides_db = filter_guides(guides_db, lenient=True) #NOTE: lenient!!
print(f"There are {len(set(guides_db.gene))} genes that have guides  out of {len(set(db_sorfs.index))} genes ({len(guides_db)} guides)")
guides_db = guides_db[guides_db.n_off_targets == 0] #so, guides in gene with no mismatches
print(f"After filtering for off-targets, There are {len(set(guides_db.gene))} genes that have guides \
     out of {len(set(db_sorfs.index))} genes ({len(guides_db)} guides)")

#guide_list = []   #can filter to ensure guides cut in first half?
#for s in db_sorfs.index:
    #guide has to be in first half, I think
#    guide_list += [(s,) + x for x in 
#                   get_all_guides(db_sorfs.chrom[s], db_sorfs.cdsStart[s], end=db_sorfs.cdsStart[s]+(db_sorfs.cdsEnd[s]-db_sorfs.cdsStart[s])//2)] 
#guides_db = pd.DataFrame(guide_list, 
#                              columns=["gene","guide_seq","chrom","cut_pos","guide_strand"])
#guides_db = filter_guides(guides_db)
#print(f"There are {len(set(guides_db.gene))} genes that have guides in the first half out of {len(set(db_sorfs.index))} genes ({len(guides_db)} guides)")

    

#actually, all the guides have GG - options for codon cut-offs are NN|N.NNN.GGN., NNN|NNN.NGG, NNN.N|NN.NNG.GNN (. = codon boundaries, | = cut)
#so we want to round to nearest boundary, up or down, and then edit second codon to TAA, and introduce a 1-bp deletion
#orrr the guide is on the opposite strand and its .CCN.NNN|NNN , NCC.NNN.N|NN , NNC.CNN.NN|N.
#so we want to round to nearest boundary, up or down, and then edit codon before previous to TAA, and introduce a 1-bp deletion
#so I think I can just change all the PAM sites specifically to - stop codon - TAA, as well as create a 1 bp deletion

#if I'm doing BOTH a stop codon and a deletion/insertion, then I need to change the code

#guides_db = guides_db[guides_db.n_off_targets == 0] #so, guides in gene with no mismatches


#NOTE: I'm changing this to include more guides
#guides targeting start
sorf_names = []
positions = []
old_windows = []
new_windows = []
guides_skipped = 0
too_early = 0
too_late = 0
cp_in_gene_list = []
tp_in_gene_list = []
same_strand_list = []
for i in guides_db.index:
    gene = guides_db.gene[i]
    gene_strand = db_sorfs.strand[gene]
    guide_strand = guides_db.guide_strand[i]
    if gene_strand == 1: #expect +/- later
        gene_strand = '+'
    else:
        gene_strand = '-'
    codons_in_gene = db_sorfs.sizes[gene]//3 #NOTE works only because no introns
    cut_pos = guides_db.cut_pos[i]
    if gene_strand == '+':
        cp_in_gene = cut_pos - db_sorfs.cdsStart[gene]
    else:
        cp_in_gene = db_sorfs.cdsEnd[gene] - cut_pos
    #actually, we'll do this differently
    #want codon to edit to be as early in sequence as possible
    #preferably intersecting PAM
    #if outside of sequence, do anything to get it inside
    #so, if guide strand is opposite of gene, target CC 
    #if guide strand is same as gene, target pos -7
    #if this lands you outside of gene, target earliest possible
    same_strand = (guide_strand == gene_strand)
    if not same_strand:
        #C[C]N NNN | NNN   ... [target]
        target_pos = cp_in_gene - 5
    else:
        # ... NN[N] NNN  NNN | NGG
        target_pos = cp_in_gene - 7
    if target_pos < 0:
        #set to 0, but confirm it can still be targeted by guide
        # cp_in_gene >= -16 not same_strand, or cp_in_gene >= -5 if same_strand
        if (cp_in_gene >= -5) or ((cp_in_gene >= -16) and (not same_strand)):
            target_pos = 0

        else:
            guides_skipped += 1
            too_early += 1
            tp_in_gene_list.append(np.nan) 
            continue

    #
    cp_in_gene_list.append(cp_in_gene)
    if target_pos >= codons_in_gene*3 - 3: #target in stop codon or beyond
        #target earliest possible for guide
        if not same_strand:
            #[C]CN NNN | NNN   ... [target]
            target_pos = cp_in_gene - 6
        else:
            #[N]N NNN NNN NNN NNN NNN | NNN NGG - target -17
            target_pos = cp_in_gene - 17
        #but now we need to confirm its in gene, set to 0 if not
        if target_pos < 0:
            target_pos = 0
        elif target_pos >= codons_in_gene*3 - 3:
            guides_skipped += 1 #even shifting it to leftmost position is not enough, give up
            tp_in_gene_list.append(np.nan) #skipped
            too_late += 1
            continue
    tp_in_gene_list.append(target_pos)

    #NOW, choose the codon that target_pos is in
    #this is floor of /3 
    codon_to_edit = int(np.floor(target_pos/3))
    assert(codon_to_edit >= 0)
    assert(codon_to_edit < codons_in_gene-1)
    sorf_names.append(gene)
    positions.append(codon_to_edit*3)
    old_windows.append(db_sorfs.seq[gene][codon_to_edit*3:codon_to_edit*3+4])
    new_window = 'TAA'
    if db_sorfs.seq[gene][codon_to_edit*3+4:codon_to_edit*3+6] == 'TG' and db_sorfs.seq[gene][codon_to_edit*3+3] != 'A': #DON'T make deletion if the deletion would create a start codon
        new_window = 'TAA'+db_sorfs.seq[gene][codon_to_edit*3+3]  #wait this only matters if 
    new_windows.append(new_window)

    
    ''' old - delete
    #assert(cp_in_gene >= 0) this isn't
    codon_to_edit = round(cp_in_gene/3) #ideally, targets PAM
    if guide_strand == '+':
        codon_to_edit += 1
    else:
        codon_to_edit += -2

    #TODO - protect against codon to edit being outside of gene?
    
    #if outside of gene [editing stop not allowed]
    #use minimum that will get it inside of gene
                                                        #this last one cannot happen because start is ATG, but still
    if codon_to_edit < 0 and (((cp_in_gene > -17) and guide_strand == '-') or (cp_in_gene >= -5)): 
        #shift to target not-pam if pam falls out of bounds
        #if guide_strand is +, already targetting farthest-right codon we can
        codon_to_edit = 0
    #what about if it falls off the other side?
    #so, theoretically this shouldn't happen [if we're targetting second half? unless super short] - min 3, 10% less than 12,
    #limit of allowed (target second-to-last)  NNN TGA NNN NNN NNN NN|N NNN GG   {but do we allow this?}
    #OR we could allow N TGA NNN NNN NNN NNN N|NN NNG G
    #lowest that would cross threshold         NNN NNN NNN NNN NNN NN|N TGA GG
    if (codon_to_edit >= codons_in_gene-1) and ((cp_in_gene <= codons_in_gene*3 + 11) and guide_strand == '+'): 
        #shift to target not-pam if pam falls out of bounds
        #if guide_strand is +, already targetting farthest-right codon we can
        codon_to_edit = max(round((cp_in_gene - 17)/3), 0) #or should I restrict this so that it is in seed??
            

    if not ((codon_to_edit < 0) or (codon_to_edit >= codons_in_gene-1)): #allow for the start codon to be edited, but not the stop
        try:
            sorf_names.append(gene)
            positions.append(codon_to_edit*3)
            old_windows.append(db_sorfs.seq[gene][codon_to_edit*3:codon_to_edit*3+4])
            new_windows.append('TAA')
        except:
            print("something went wrong")
            print(guides_db.loc[i])
            print(f"gene has {codons_in_gene} codons, cp_in_gene is {cp_in_gene}, codon_to_edit is {codon_to_edit}")
            raise ValueError()
    else:
        guides_skipped += 1
        '''

guides_db["cp_in_gene"] = cp_in_gene_list
guides_db["tp_in_gene"] = tp_in_gene_list
guides_db.to_csv('tmp/guides_db.csv')


print(f"Had to skip {guides_skipped} guides becuase edit would have been outside of gene, {too_early} of which were upstream of start, {too_late} downstream of end")
print(f"we have {len(positions)} target positions covering {len(set(sorf_names))} genes")

#only the ... first 3? 
#variant_inserts = make_library(sorf_names, positions, old_windows, deletion=True, w=2)

variant_inserts = make_library(sorf_names, positions, old_windows, new_windows, w=4, seed_size=17, lenient=True)
variant_inserts.to_csv('tmp/unfiltered_inserts.csv')

print('Length variants:', len(variant_inserts))
#print('Variants targeting 1st position:', np.sum(variant_inserts.w_gene_pos == 0))
#print('Highest position targeted (in codons)', max(variant_inserts.w_gene_pos)//3)
print('Unique genes targeted:', len(set(variant_inserts.gene)))

variant_inserts['gene_length'] = [db_sorfs.cdsEnd[i] - db_sorfs.cdsStart[i] for i in variant_inserts.gene]
x = len(set(variant_inserts[((variant_inserts.w_gene_pos / variant_inserts.gene_length) <= 0.5)].gene))
print('Genes with edit position less than halfway down gene', x)
x = len(set(variant_inserts[((variant_inserts.w_gene_pos / variant_inserts.gene_length) <= 0.25)].gene))
print('Genes with edit less than 25 percent down the gene', x)
x = len(set(variant_inserts[(variant_inserts.w_gene_pos == 0)].gene))
print('Genes with the first start codon targeted', x)

variant_inserts['percent_loc'] = [percent_loc(variant_inserts.gene[i], variant_inserts.w_codon_num[i]) for i in variant_inserts.index]
variant_inserts['downstream_starts'] = [downstream_starts(variant_inserts.gene[i], variant_inserts.w_codon_num[i]) for i in variant_inserts.index]
variant_inserts.to_csv('brar_data/sorfs_1_lenient.csv')

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

