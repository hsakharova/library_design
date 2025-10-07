import pandas as pd
import numpy as np
import re
from ast import literal_eval
import collections
import os
print('restarted')

#load_files
db_genes = pd.read_csv("data/db_genes.csv", sep=",", index_col=0, header=0)
db_genes = db_genes.where(pd.notnull(db_genes), None) 
#print(db_genes.exonStarts.dtype, type(db_genes.exonStarts[0]), db_genes.exonStarts[0])
db_genes.exonStarts = [literal_eval(x) for x in db_genes.exonStarts]
db_genes.exonEnds = [literal_eval(x) for x in db_genes.exonEnds]
codon_info = pd.read_pickle("data/codon_info.pkl")
index_file = 'yeast_genes/sacCer_index'


#TODO even more hacky
reference_sacCer3 = True

#TODO hacky
def set_db_genes(db):
    global db_genes 
    db_genes = db
    print(f'Set db_genes of {len(db_genes)} genes, first gene is {db_genes.index[0]}')

def get_db_genes():
    return db_genes

def set_index_file(fname):
    global index_file
    index_file = fname
    print(f'Set index as {fname}')


#TODO: unhack this
def get_chrom_dict(fname='brar_data/sk1_original.fasta'):
    chrom_dict  = {}
    with open(fname, 'r') as f:
        lines = f.readlines()
        chrom = "no_chrom"
        for i, line in enumerate(lines):
            line =line.strip()
            #if i <= 10:
                #print(line)
            if line[0] == '>':
                chrom = line[1:]
                assert(not chrom in chrom_dict)
                chrom_dict[chrom] = ""
            else:
                chrom_dict[chrom] += line
        print(f'Using {fname} as genome. Chromosome lengths:')
        print([(k, len(chrom_dict[k])) for k in chrom_dict])
    return chrom_dict

#chrom_dict = get_chrom_dict()

def set_chrom_dict(new_chrom_dict):
    global chrom_dict
    chrom_dict = new_chrom_dict
    global reference_sacCer3
    reference_sacCer3 = False

def obtain_chrom_dict():
    print('reference is', reference_sacCer3)
    if not reference_sacCer3:
        return chrom_dict

def get_seq(chrom, start, end, original=True):
    if not reference_sacCer3:
        return chrom_dict[chrom][start:end]
    else:
        with open("yeast_genes/" + chrom + ".fa", "r") as f:
            lines = f.readlines()
            seq = "".join(lines[1:])
            seq = "".join(seq.split("\n"))
        return seq[start:end]



'''
User-facing function. Generates a CRISPEY library of donor-guide inserts targeting specified windows in the given genes.
If deletion is True, donor-guide inserts will delete specified windows. Otherwise, specified windows will be converted into
synonymous sequences using the fastest-possible codons.

Inputs:
- genes: a list of strings, the yeast ORF names of the genes being targeted
- positions: a list of integers, the positions (in base pairs) of the targeted windows inside each gene (0-indexed)
- original_windows: a list of strings, the original DNA nucleotide sequence of the windows being targeted. 
    (Note that genes, positions, and original_windows should be the same length)
- deletion: boolean. Whether or not the generated CRISPEY guide-donor inserts should delete the target windows (True)
    or replace the target windows with synonymous sequences of the fastest-possible synonymous codons. (False)
- w: size of windows (in bp). If deletion is False, must be a multiple of 3 (to target codons)
- do_not_filter: boolean. Whether or not to filter the generated guides. The filtering steps (by default) are to remove all 
    guides with off-target matches elsewhere in the yeast genome (guides that have an off-target match with less than 3 mismatches).
    Off-target matches are found using bowtie2 (installation of bowtie2 is required). Donor-guide inserts are then selected to 
    not have long repeats and to have an edited nucleotide within the seed region of the guide. 
 TODO added new_windows
Outputs:
- inserts_db: Pandas dataframe containing the donor-guide CRISPEY insert sequences along with relevant information 
'''
def make_library(genes, positions, original_windows, new_windows=None, w=3, do_not_filter=False, seed_size=7, lenient=False, prefix='v_'):
    assert(len(genes)==len(positions))
    if new_windows is not None:
        assert(len(set(zip(genes, positions, new_windows))) == len(positions)), 'targets should be unique' #
    else:
        assert(len(set(zip(genes, positions))) == len(positions)), 'targets should be unique'
    print("Number of targets:", len(genes))
    idb_list = []
    for i in range(len(genes)):
        #give an update every 10% of the way
        if len(genes) < 10 or i%(len(genes)//10)==0:
            print("processing target", i, "out of", len(genes))
        gene = genes[i]
        pos = positions[i]
        #if convert:
        #    pos = convert_pos_to_cn(pos, gene, w=w) 
        if gene_pos_to_chrom_pos(gene,pos,w=w, throw_error=False, w_in_bp=True) is None:
            #skip if window not fully in one exon,
            #make inserts would throw error.
            continue
        if new_windows is None:
            new_window = None
        else:
            new_window = new_windows[i]
        idb = make_inserts(gene, pos, original_windows[i], new_window = new_window, w=w)
        #sanity_checks(idb, w=w) for speed, will do later
        # add stdev info, too.
        idb_list.append(idb)
    inserts_db = pd.concat(idb_list)
    inserts_db = inserts_db.reset_index()
    inserts_db = inserts_db.drop("index", axis=1) 
    #CHANGED: add new_window to target
    inserts_db["target"] = [inserts_db.gene[i]+"_"+str(inserts_db.w_gene_pos[i]) + '_' + inserts_db.new_window[i] 
                          for i in inserts_db.index]
    #CHANEGD: add new_window to name
    inserts_db["name"] = [prefix+inserts_db.gene[i]+"_"+str(inserts_db.w_gene_pos[i]) + '_' + inserts_db.new_window[i] +'_guide_'+
                          inserts_db.chrom[i]+'_'+inserts_db.guide_strand[i]+'_'+str(inserts_db.cut_pos[i])
                          for i in inserts_db.index]
    inserts_db = inserts_db.set_index('name')
    original_length = len(inserts_db)
    #drop duplicates - TODO - how would duplicates happen in the first place? unsure but they did - investigate
    inserts_db = inserts_db.drop_duplicates()
    print(f"Dropped {original_length - len(inserts_db)} duplicates")
    if do_not_filter:
        print("DID NOT FILTER")
        return inserts_db
    #Ensure no 'N' in insert seq
    inserts_db = inserts_db[[True if 'N' not in x else False for x in inserts_db.insert_seq]]
    print(f"After removing any guide-donors with 'N': {len(inserts_db)} guides \
          for {len(set(list(inserts_db.target)))} targets") 
    #instead of dropping mismatches/etc., lets keep them?
    #filter for change
    inserts_db = filter_for_change(inserts_db)
    print(f"After removing any guide-donors that would target themselves: {len(inserts_db)} guides \
          for {len(set(list(inserts_db.target)))} targets") 
    print("")
    print("Before filtering for off-target matches:", len(inserts_db), 
          "guides for", len(set(list(inserts_db.target))), "targets in",
          len(set(inserts_db.gene)),"genes")
    print("Filtering guides, allow offtargets with any mismatch in seed", lenient)
    inserts_db = filter_guides(inserts_db, lenient=lenient) 
    print("Off targets:", sorted([x for x in collections.Counter(inserts_db.n_off_targets).items() if x[0]<=10]))
    print("Guides with more than 10 off targets:", np.sum(inserts_db.n_off_targets > 10))
    inserts_db = inserts_db[inserts_db.n_off_targets == 0]
    print("")
    print("After filtering for off-target matches:", len(inserts_db), 
          "guides for", len(set(list(inserts_db.target))), "targets in",
          len(set(inserts_db.gene)),"genes")
    print("running sanity checks")
    sanity_checks(inserts_db, additional_checks = True, synonymous=(new_windows is None))
     
    add_edits_info(inserts_db)
    inserts_db = count_replicates(inserts_db)
    select_guides(inserts_db, seed_size=seed_size)
    #TODO: do I want to ensure that guides are unique? nah
    inserts_db = inserts_db[inserts_db.selected]
    print("")
    print("After selecting for repeats/edits in guide:", len(inserts_db), 
          "guides for", len(set(list(inserts_db.target))), "targets in",
          len(set(inserts_db.gene)),"genes")
    inserts_db = count_replicates(inserts_db)
    inserts_db = inserts_db.reset_index()
    return inserts_db #filtered for selected 

''' 
Make all possible CRISPEY oligomer inserts targeting a window of length w starting 
at a position (gene_pos) within a gene. If deletion is True, the oligomer inserts will
create a deletion of length w at position gene_pos within the gene. Otherwise, the inserts
will replace all slow codons within the window with the fastest possible codons. 

h_arm_length sets the length of the homology arms, and original_window provides the original
sequence of the window to be modified (primarily for downstream sanity tests)
'''
def make_inserts(gene, gene_pos, original_window, new_window=None, w=3, h_arm_length=50, alt_fast=True):
    #TODO MAKING CHANGES
        #gene_pos, w, and h_arm_length are in nucleotides
        #if deletion is False, change a window of w/3 codons to the fastest possible codons
        #if deletion is True, delete a window of w bp
    #if None, replace  window of slow codons with fast codons
    #else, replace original_window with new_window
    #so, no deletion
    
    if new_window is None:
        if gene_pos % 3 != 0:
            raise ValueError('Gene pos is not in frame')
        if w % 3 != 0:
            raise ValueError('Window is not mod 3 - out of frame')
    else:
        assert(len(original_window) >= len(new_window)) #NO INSERTIONS ALLOWED FOR NOW
    #    if new_window is not None:
    #        assert((len(new_window)==len(original_window))), 'New and old windows must be same length'
    #else:
    #    assert((new_window == '') or (new_window == None)), 'changing window unsupported for deletion'
    assert(len(original_window) == w) #w is the length of the original window
    constant_5ph = "GAGTTACTGTCTGTTTTCCT"
    constant_rt = "AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCA"
    constant_3ph = "GTTTCAGAGCTATGCTGGAA"
    chrom, chrom_pos, gene_strand = gene_pos_to_chrom_pos(gene,gene_pos,w=w, w_in_bp=True) #TODO: why do you need w here?
    guide_list = get_all_guides(chrom, chrom_pos, w=w) #w for get all guides is in bp
    inserts_db = pd.DataFrame(guide_list, 
                              columns=["guide_seq","chrom","cut_pos","guide_strand"])
    inserts_db["gene"] = gene
    inserts_db["gene_strand"] = gene_strand
    inserts_db["w_codon_num"] = gene_pos//3 #in codons (round down)
    inserts_db["w_gene_pos"] = gene_pos #base pairs into gene.
    #start of window as position in gene. does NOT include introns [cds sequence]
    inserts_db["w_chrom_pos"] = chrom_pos #position in chromosome
    inserts_db["w_size"] = w #in bp, not codons!
    inserts_db["target_window"] = original_window
    #I would like this all to be in a pandas df 
    #with gene, gene_pos, chrom_pos, strand, 
    old_window = get_seq(chrom, chrom_pos, chrom_pos+w) 
    if gene_strand == "-":
        old_window = rev_complement(old_window)
    inserts_db["old_window"] =  old_window #this should be the same as target window - keep both as a sanity check
    if new_window is None:
        #new_window = "" #if deletion - deletion now specififed with ""
        #if not deletion: #change to a new window
        new_window = slow_to_fast(old_window) #change all slow codons to fast codons
        #OR use this alternative set of fast codons
        #mostly because these are almost the same er
        #but higher cf/tai
        #note that as written this will only change single codons, not a window of multiple codons
        if alt_fast:
            slow_to_alt_fast_dict = {"CGG":"AGA", "CGA":"AGA", "AGG":"AGA", "CGC":"AGA", #ARG
                                    "CCG":"CCA", "CCC":"CCA", #PRO
                                    "GCG":"GCT", "GCA":"GCT"} #ALA
            if old_window in slow_to_alt_fast_dict.keys():
                new_window = slow_to_alt_fast_dict[old_window]
    #else we are keeping the given new_window
    inserts_db["new_window"] = new_window
    wt_donor_list = []
    mod_donor_list = []
    insert_list = []
    change = len(original_window) - len(new_window)
    for i in inserts_db.index:
        #donor direction --> should match guide direction ??
        cp = inserts_db.cut_pos[i]
        guide_strand = inserts_db.guide_strand[i]
        guide_seq = inserts_db.guide_seq[i]
        #if not deletion:
            #if no indels, just substitution
        #    donor = get_seq(chrom, cp-h_arm_length, cp+h_arm_length)
        #    loc_w = h_arm_length + chrom_pos - cp
        #if deletion:
            #depending on if chrom_pos is before/after cut_pos
        shift = max(min(cp - chrom_pos, change), 0) #so between del_size (w) and 0
        donor = get_seq(chrom, cp-h_arm_length-shift, cp+h_arm_length+change-shift) #+ strand
        loc_w = h_arm_length + chrom_pos - cp + shift   #TODO what is loc_w?
        window_to_insert = new_window
        if gene_strand == "-":
            window_to_insert = rev_complement(window_to_insert)
        mod_donor = donor[:loc_w]+window_to_insert+donor[loc_w+w:]  #+ strand
        if guide_strand == "-":
            donor=rev_complement(donor)
            mod_donor = rev_complement(mod_donor) #So if both are -, window will be normal again
        wt_donor_list.append(donor)
        mod_donor_list.append(mod_donor)
        insert_list.append(constant_5ph + mod_donor + constant_rt 
                           + guide_seq + constant_3ph)
    inserts_db["wt_seq"] = wt_donor_list
    inserts_db["donor"] = mod_donor_list
    inserts_db["insert_seq"] = insert_list
    #NOTE: SO THIS MEANS THAT DONOR IS IN THE DIRECTION OF THE GUIDE
    #NOT necessarily in the direction of the gene
    return inserts_db


#so I guess now, if we want all [?] possible ways to knock down a gene
#or do we want to filter some preemptiveley - e.g. early in a sequence? multiple guides targetting same site?
#want to introduce stop+frameshift ...?
#maybe, first, check which genes have guides that work




'''Returns all possible guide sequences that will intersect with a given window
Window can be defined by start:start+w (note that start is with respect to chromosome, not strand)
or by start:end'''
def get_all_guides(chrom,start,end=None, w=3): #20 is already too far away, actually - but whatever
    #w and max_dist are in bp
    if end is None:
        end = start+w
    guide_length = 20
    #back_dist = max_dist-w-4 #distance to search back 
    #(4-> 'cut' 1 2 3 N | G G --> to limit where to look for GG)
    back_dist = 1 #only | G  w_0/G w_1 w_2 .... is acceptable - where to look for GG
    #fw_dist = max_dist-w+6 #distance to search fw
    #(4-> 'cut' 1 2 3 N G G | --> to limit where to look for GG)
    fw_dist = guide_length + 2 #(fw_dist is from end) (-1 from end so that last w intersects, +23 guide including PAM = 22)
    #so  w_0 w_1 |w_2/-20 -19 ....... -1 N G G| is acceptable
    #setting so that if window intersects at all, it's ok
    fw_seq = get_seq(chrom, start-back_dist, end+fw_dist)
    rev_seq = rev_complement(get_seq(chrom, start-fw_dist, end+back_dist))
    fw_indices = [m.start() + start - back_dist - 4 #to get location of cut
                  for m in re.finditer('GG',fw_seq)]
    rev_indices = [end -m.start() + back_dist + 4 #to get location of cut
                   for m in re.finditer('GG', rev_seq)]
    #print(rev_seq)
    #print([m.start() for m in re.finditer('GG', rev_seq)])
    fw_guides = [(get_guide(chrom,cp,"+"),chrom,cp,"+") 
                 for cp in fw_indices]
    rev_guides = [(get_guide(chrom,cp,"-"),chrom,cp,"-")
                 for cp in rev_indices]
    return fw_guides + rev_guides




"Function to check that database of inserts (idb) was constructred correctly"
def sanity_checks(idb, additional_checks=False, synonymous=False): #assumes 50 bp arms, cut pos in middle
    for i in idb.index:
        w = idb.w_size[i] #IN BP
        assert(len(idb.donor[i])==100)
        shift = 0
        del_size = len(idb.old_window[i]) - len(idb.new_window[i])
        guide_s = 1
        if idb.guide_strand[i]=="-":
            guide_s = -1
        shift = min(max(idb.cut_pos[i]-idb.w_chrom_pos[i], 0), del_size)
        if idb.guide_strand[i]=="-":
            shift = del_size-shift #I think?
        #assert(len(idb.wt_seq[i])+del_size==100)
        assert(len(idb.wt_seq[i])-del_size==100)
        assert(idb.wt_seq[i][33+shift:53+shift] == idb.guide_seq[i]), str(idb.gene[i]) + "_" + str(idb.w_codon_num[i])
        assert(idb.wt_seq[i][54+shift:56+shift]=="GG") #PAM assertion 
        assert(idb.target_window[i]==idb.old_window[i]), i
        assert(len(idb.insert_seq[i])==194)
        '''
        if not deletion:
            w_on_strand = idb.new_window[i]
            donor_on_strand = idb.donor[i]
            if idb.gene_strand[i] == "-":
                w_on_strand = rev_complement(w_on_strand)
            if idb.guide_strand[i] == "-":
                donor_on_strand = rev_complement(donor_on_strand)
            assert(idb.guide_seq[i] in idb.wt_seq[i])
            #assert(idb.wt_seq[i].index(idb.guide_seq[i]) == 33) #this could fail if guide in sequence twice??
            assert(w_on_strand in donor_on_strand)
            #assert(donor_on_strand.index(w_on_strand) == 50 + idb.w_chrom_pos[i] - idb.cut_pos[i]) #this doesn't work if multiple same windows, eg w=1
            try:
                assert(donor_on_strand[50+idb.w_chrom_pos[i]-idb.cut_pos[i]:50+idb.w_chrom_pos[i]-idb.cut_pos[i]+w]==w_on_strand)
            except:
                print("Something went wrong")
                print(idb.gene[i])
                print(idb.gene_strand[i])
                print(idb.w_chrom_pos[i], idb.cut_pos[i], w, idb.w_gene_pos[i])
                print(donor_on_strand)
                print(w_on_strand)
                print(donor_on_strand[50+idb.w_chrom_pos[i]-idb.cut_pos[i]:50+idb.w_chrom_pos[i]-idb.cut_pos[i]+w])
                raise ValueError()
            nm_window = num_mismatch(idb.old_window[i], idb.new_window[i])
            nm_donor = num_mismatch(idb.wt_seq[i], idb.donor[i])  #TODO: target window *is* old window, not 
            assert(nm_window==nm_donor)
            assert(len(idb.new_window[i])==w)
        else:
        '''
        insert_pos = 50 + (idb.w_chrom_pos[i] - idb.cut_pos[i])*guide_s + shift
        if guide_s == -1:
            insert_pos = insert_pos - del_size - len(idb.new_window[i]) #TODO: correct?
        old_seq = idb.old_window[i]
        assert(len(old_seq)==del_size+len(idb.new_window[i]))
        if idb.guide_strand[i] != idb.gene_strand[i]:
            old_seq = rev_complement(old_seq)
        assert_msg = idb.gene[i] + idb.gene_strand[i] + idb.guide_strand[i] + " "
        assert_msg += old_seq + " " + str(shift) +  " " + str(insert_pos) + " "+ idb.donor[i] + " " + idb.wt_seq[i]
        assert(idb.donor[i][:insert_pos]+old_seq+idb.donor[i][insert_pos+len(idb.new_window[i]):]==idb.wt_seq[i]), assert_msg
        assert(idb.wt_seq[i][insert_pos:insert_pos+del_size+len(idb.new_window[i])]==old_seq), assert_msg
        if additional_checks:
            msg = [idb.gene[i], idb.w_codon_num[i], idb.w_size[i]]
            #check that window does not cross introns/ends
            assert(not intersecting_intron(idb.gene[i], idb.w_chrom_pos[i], idb.w_chrom_pos[i]+w)), msg
            #check that slow is indeed correct
            assert(get_seq_in_gene(idb.gene[i], idb.w_gene_pos[i], idb.w_gene_pos[i]+w) == idb.target_window[i])
            if idb.gene_strand[i] == '-':
                assert(rev_complement(get_seq(idb.chrom[i], idb.w_chrom_pos[i], idb.w_chrom_pos[i]+w)) == idb.target_window[i])
            else:
                assert(get_seq(idb.chrom[i], idb.w_chrom_pos[i], idb.w_chrom_pos[i]+w) == idb.target_window[i])  
            #a change was made!
            assert(idb.old_window[i] != idb.new_window[i])
            assert(idb.wt_seq[i] != idb.donor[i])
            #TODO - check guide seq + NGG not in donor --
            assert(not re.search(idb.guide_seq[i] + '.GG', idb.donor[i]))

            #check that amino acid sequence slow, fast is the same
            #if not deletion:
            if synonymous:
                assert(amino_acid_seq(idb.old_window[i]) == amino_acid_seq(idb.new_window[i]))
                #fast is faster than slow?  #ah but for neutral [1 codon change] it's not actually faster .... oops
                if idb.w_size[i] != 3: #the following is not true for single codon changes (where I used alt_fast)
                    assert(slow_to_fast(idb.old_window[i]) == idb.new_window[i]) #
            assert(idb.w_codon_num[i] == idb.w_gene_pos[i] // 3) #I think this is correct
            #TODO slow/fast are misnomers for neutral
            #check that constant sequences are correct, and donor/guide in correct spot
            assert(idb.insert_seq[i][:20] == "GAGTTACTGTCTGTTTTCCT")
            assert(idb.insert_seq[i][-20:] == "GTTTCAGAGCTATGCTGGAA")
            assert(idb.insert_seq[i][-74:-40] == "AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCA")
            assert(idb.insert_seq[i][-40:-20] == idb.guide_seq[i])
            assert(idb.insert_seq[i][20:120] == idb.donor[i])
            #check that from w_codon_num get correct w_chrom_pos
            chr_ps = get_chrom_pos(idb.gene[i],idb.w_gene_pos[i],idb.w_size[i])
            assert(chr_ps == idb.w_chrom_pos[i]), [chr_ps, idb.w_chrom_pos[i], idb.gene[i], 
                                                   idb.w_gene_pos[i], idb.w_size[i]]

'''Function to recover sequences from fasta files corresponding to yeast chromosomes
- substitute if sequences must be recovered from elsewhere'''   
'''      
def get_seq(chrom,start,end):
    with open("yeast_genes/" + chrom + ".fa", "r") as f:
        lines = f.readlines()
        seq = "".join(lines[1:])
        seq = "".join(seq.split("\n"))
    return seq[start:end]
'''

'''
Helper functions below
'''
def rev_complement(seq):
    nt_pairs = {"A":"T","G":"C","C":"G","T":"A","N":"N","-":"-"} #do i want "-" to be here???
    rev_seq = seq[::-1]
    rev_seq = [nt_pairs[nt] for nt in rev_seq]
    return "".join(rev_seq)

def get_guide(chrom, cut_pos, strand):
    #TODO --> assert NGG??
    if strand=="-":
        return rev_complement(get_seq(chrom,cut_pos-3,cut_pos+17))
    else:
        return get_seq(chrom, cut_pos-17, cut_pos+3)
    
def guide_list_to_dict(gl):
    return {str(x[1])+"_"+str(x[2])+x[3]:x[0] for x in gl}

#guides_dict is guide_name : guide_seq
def make_fastq_file(fname, guides_dict, add_pam=False, ignore_bases = 0):
    #if add_pam, we want to add "NGG" for the guide
    #print(guides_dict)
    pam = ""
    if add_pam:
        pam ="NGG" #question --> note that this 'N' counts as an extra mismatch
    with open(fname, 'w') as f:
        for key in guides_dict:
            try:
                f.write(">"+key+"\n")
                f.write(guides_dict[key][ignore_bases:]+pam+"\n")
            except:
                print(f"make_fastq_file failed to write")
                print(f"key is {key} and guides_dict[key] is {guides_dict[key]}")
                raise ValueError()
    return

#position in gene to position in chromosome
## NOGAP pos in gene
#TODO - test
def gene_pos_to_chrom_pos(gene, pos, w=10, throw_error=True, w_in_bp=False):
    #remember, pos is start of window relative to gene
    #but chrom_pos needs to be _start_ of window on + strand
    #(so if gene is on - starnd its actually the end of the window)
    #pos are in bp, here
    if not w_in_bp:
        w = w*3
    chrom = db_genes.chrom[gene]
    strand=db_genes.strand[gene]
    if strand == 1:
        strand = "+"
    elif strand == -1:
        strand = "-"
    num_exons = db_genes.exonCount[gene]
    #TODO: pretty sure start/end in db_genes refers to
    #loc on chrom (not strand) but check
    exon_starts = db_genes.exonStarts[gene]
    exon_ends = db_genes.exonEnds[gene]
    if strand=="+":
        gene_start = db_genes.cdsStart[gene]
        chrom_pos = gene_start + pos
        if num_exons != 1:
            intron_bp = 0
            intron_lengths = [exon_starts[i+1] - exon_ends[i]
                              for i in range(num_exons-1)]
            for i in range(num_exons):
                if chrom_pos + intron_bp >= exon_ends[i]:
                    intron_bp += intron_lengths[i]
                else:
                    #ensure window fully in exon.
                    assertion_error = str(chrom_pos)+" "+str(intron_bp)+" "+str(exon_starts[-i-1])
                    assertion_error+= " "+str(gene)+" "+str(pos)+" "+strand
                    assertion = (chrom_pos + intron_bp + w < exon_ends[i])
                    if throw_error:
                        assert(assertion), assertion_error
                    elif not assertion: #if not throwing error, return None
                        print("Not fully in exon "+assertion_error)
                        #it's true! sometimes windows aren't fully in exons.
                        #can't use these
                        #will return None
                        return None
            chrom_pos = chrom_pos + intron_bp
    if strand=="-": #TODO - check
        gene_start = db_genes.cdsEnd[gene]
        chrom_pos = gene_start - pos #on + str it is end of window
        if num_exons != 1:
            intron_bp = 0
            intron_lengths = [exon_starts[i+1] - exon_ends[i]
                              for i in range(num_exons-1)][::-1] #reverse!
            for i in range(num_exons):
                if chrom_pos - intron_bp <= exon_starts[-i-1]:
                    intron_bp += intron_lengths[i]
                else:
                    #ensure window fully in exon
                    assertion_error = str(chrom_pos)+" "+str(intron_bp)+" "+str(exon_starts[-i-1])
                    assertion_error+= " "+str(gene)+" "+str(pos)+" "+strand
                    assertion = (chrom_pos - intron_bp - w >= exon_starts[-i-1])
                    if throw_error:
                        assert(assertion), assertion_error
                    elif not assertion:
                        print("Not fully in exon "+assertion_error)
                        return None
            chrom_pos = chrom_pos - intron_bp
        chrom_pos = chrom_pos - w #NOW it is start of window on + strand
    return chrom, chrom_pos, strand
    #start should be start of 10 bp window on + strand 

def num_mismatch(s1,s2):
    assert(len(s1)==len(s2))
    return sum([1 for i in range(len(s1)) if s1[i] != s2[i]])

#TODO: which one?
def get_chrom_pos(gene, gene_pos, w): #another version exists, but I wanted this for test. w is in bp
    exon_starts = db_genes.exonStarts[gene]
    exon_ends = db_genes.exonEnds[gene]
    exon_lengths = list(np.array(exon_ends) - np.array(exon_starts))
    strand = db_genes.strand[gene]
    if strand == -1:
        exon_lengths = exon_lengths[::-1]
        exon_starts_ = exon_ends[::-1]
        exon_ends = exon_starts[::-1]
        exon_starts = exon_starts_
    to_subtract = [0] + exon_lengths 
    assert(gene_pos >= 0)
    assert(len(to_subtract) >= 1)
    assert(to_subtract[0] == 0)
    assert(gene_pos <= sum(exon_lengths))
    pos_after_exon_start = [gene_pos - sum(to_subtract[:i]) for i in range(len(to_subtract)) 
                            if gene_pos-sum(to_subtract[:i]) >=0] #cut off
    j = len(pos_after_exon_start) - 2
    assert(j >= 0)
    assert(j < len(exon_starts)), str(exon_starts) + str(exon_lengths)+gene+"_"+str(gene_pos)+"_"+str(j)
    if strand == -1:
        return exon_starts[j]-pos_after_exon_start[-1]-w
    else:
        return exon_starts[j]+pos_after_exon_start[-1]
    
def amino_acid_seq(seq):
    assert(len(seq)%3 == 0)
    return [codon_info.AA[seq[i*3:i*3+3]] for i in range(len(seq)//3)]
        
def get_seq_in_gene(gene, start_bp, end_bp):
    assert(start_bp >= 0)
    assert(start_bp < end_bp)
    exon_starts = db_genes.exonStarts[gene]
    exon_ends = db_genes.exonEnds[gene]
    assert(len(exon_starts)==len(exon_ends))
    chrom = db_genes.chrom[gene]
    seq = "".join([get_seq(chrom, exon_starts[i], exon_ends[i]) for i in range(len(exon_starts))])
    if db_genes.strand[gene] == -1:
        seq = rev_complement(seq)
    assert(end_bp <= len(seq))
    return seq[start_bp:end_bp]
        
        
def round_abs(x):
    if x >= 0:
        return np.ceil(x)
    return np.floor(x)

def intersecting_intron(gene, start, end, verbose=False):
    #TODO - also make sure it's not interesecting start/stop of gene!!!
    #start and end refer to chrom pos on + strand
    assert(start < end)
    exon_starts = np.array(db_genes.exonStarts[gene])
    exon_ends = np.array(db_genes.exonEnds[gene])
    if start < exon_starts[0] or end > exon_ends[-1]:
        if verbose:
            print("out of bounds")
        return True  #out of bounds
    if np.where(exon_ends >= end)[0][0] in np.where(exon_starts <= start)[0]:
        if verbose:
            print(np.where(exon_starts<=start)[0], np.where(exon_ends >= end)[0], "not intersect")
        return False
    else:
        if verbose:
            print(np.where(exon_starts<=start)[0], np.where(exon_ends >= end)[0], "intersect")
            print(start, end, exon_starts, exon_ends)
        return True

def filter_for_change(inserts_db):
    #this HAS to be used if mutations are just deletions without accounting for context
    #checks that the guide_seq + NGG has NO match in the donor
    changed = np.array([not bool(re.search(inserts_db.guide_seq[i] + '.GG', inserts_db.donor[i])) for i in inserts_db.index]).astype(bool)
    return inserts_db[changed]

#lenient allows for any mismatch in seed
def filter_guides(inserts_db, guide_col="guide_seq", lenient=False, random_guides=False):
    #filter guides for those that have no off-targets (<= 3 mismatches) 
    #add 'n_off_targets' to 'filter_guides'
    save_f = "tmp/inserts.fa"
    #TODO
    #seed here will be considered the 10 nt upstream of the guide
    make_fastq_file(save_f, {str(i):inserts_db[guide_col][i] for i in inserts_db.index}, add_pam=False, ignore_bases=0)
    #%cd /mnt/lareaulab/shelen/slow_codons_2021/data/
    #os.system(f"bowtie -f -v 3 -a {index_file} tmp/inserts.fa > tmp/bowtie_out.csv")
    os.system(f"bowtie -f -v 3 -a {index_file} tmp/inserts.fa > tmp/bowtie_out.csv") 
    #wanted to change to -v 4 because 'N' in PAM is always a mismatch
    #unfortunate. Turns out you can't have mismatch-4, only up to 3
    #I coulddd align with PAM, 3 (actually 2) mismatches, drop all where there is no mismatch in seed or PAM?
    #this would still "allow" any with 3 mismatches in not-seed-or-pam 
    #should I allow this...?
    #I could be harsher, and just align the seed and PAM, with v-2 .... this will remove guides with tons of mismatches outside the seed
    #but allow them as long as there is at least 1 difference in seed or PAM
    #The only other solution I can think of is somehow combining info from several bowtie runs?
    #Ok FOR NOW we will allow guides that have .... any mismatches in the seed region or PAM? so v-1
    #to be slightly stricter we can allow any mismatches in PAM, or 2 mismatches in seed (10)
    #could try any mismatches in PAM, or 1 mismatch in seed (7? 8?)
    #%cd /mnt/lareaulab/shelen/slow_codons_2021/code/
    #naahhhh these are too short. 
    #I think you can first do no-PAM, up to 3 mismatches
    #... and then check all the 3-mismatches for if they have a PAM sequence after them...?
    #... or do 
    load_f = "tmp/bowtie_out.csv"
    bowtie_info = pd.read_csv(load_f, sep="\t", header=None, names=["guide_id","strand","chr","chr_pos","seq","aln","?","mismatch"])
    #TODO - allow mismatches if the mismatch is in the PAM region? [21,22]
    #bowtie_info = bowtie_info[~(bowtie_info.mismatch.str.contains('21') | bowtie_info.mismatch.str.contains('22'))]
    bowtie_info = bowtie_info.fillna("") #some mismatch entries are empty
    if lenient:
        #if lenient, treat any with a mismatch in the seed as a mismatch [seed region = 7 bp]
        searchfor =  '|'.join([str(x) for x in range(13,20)])
        print(f"Removing {np.sum(bowtie_info.mismatch.str.contains(searchfor))} off-targets with mismatch in seed")
        #print(bowtie_info.mismatch.str.contains(searchfor)[0:10])
        #print(bowtie_info.mismatch.str.contains(searchfor).dtype)
        bowtie_info = bowtie_info[~(bowtie_info.mismatch.str.contains(searchfor))]
    #TODO - allow mismatches if there are off-targets in the seed region?
    #TODO - check if there is a PAM after the sequence?
    print('Checking for PAMs')
    #need to account for strand       
    bowtie_info["PAM"] = [((bowtie_info.strand[i]=='+') and (get_seq(bowtie_info.chr[i], bowtie_info.chr_pos[i]+21, bowtie_info.chr_pos[i]+23) == 'GG')) or
                          ((bowtie_info.strand[i]=='-') and (get_seq(bowtie_info.chr[i], bowtie_info.chr_pos[i]-3, bowtie_info.chr_pos[i]-1) == 'CC')) 
                           for i in bowtie_info.index]
    if not random_guides:
        bowtie_info = bowtie_info[bowtie_info.PAM] #only those that have a PAM are valid alignments 
        #for random guides, let's get guides that align nowhere at all ... ?
    #maybe there's a lot of guides that only align to one spot, but have multiple targets
    #TODO: this works now, clean it up
    inserts_db["n_off_targets"] = [sum(bowtie_info.guide_id == i)-1 for i in inserts_db.index]
    if not random_guides:
        print(inserts_db[inserts_db.n_off_targets < 0])
        assert(np.all(inserts_db.n_off_targets >= 0)), "some guides not aligned by bowtie at all??"
    return inserts_db
    #guides_to_use = [seq for seq in inserts_db[guide_col] if sum(bowtie_info.seq==seq) == 1]
    #return guides_to_use #sequences of guides

def add_edits_info(idb): 

    #assumes 50 bp arms
    #TODO -- this is hacky in how it handles deletion, fix.
    #how was this handled before?
    
    #add info about edits to db
    #how many edits away from wt
    #how many edits in guide
    #position of edit closest to cut site?
    # -17 .... -1 |cut| 0 1 2   (indexing scheme)
    #change indexing scheme to:  -4 |cut| -3 -2 -1 0 1 2 (where 1 and 2 are PAM GG)
    #longest homopolymer in [mutated] donor 
    #TODO - so -- how does this work with direction of guide and donor?
    #TODO - does this work for deletions....?? or do you need to shift...?
    n_edits = []
    n_edits_guide = []
    edits_in_pam = []
    edit_pos = []
    homopolymer = []
    for i in idb.index:
        #this may overemphasize impact of deletion - basically says "oh if there was a deletion and everything shifted, everything is a mismatch"
        #calculate where PAM should be now given deletion [if it was the PAM being deleted, consider everything to have 'shifted' over]
        #and then based on that PAM, count number of mismatches.

        #OR! we could do both - also say which were deleted?
        shift = 0
        del_size = len(idb.old_window[i]) - len(idb.new_window[i])
        shift = min(max(idb.cut_pos[i]-idb.w_chrom_pos[i], 0), del_size)
        if idb.guide_strand[i]=="-":
            shift = del_size-shift #I think?
        #so, in donor, centered around cut_pos
        #cut_pos is 50
        #TODO - is this assumption valid?
        wt = idb.wt_seq[i]
        donor = idb.donor[i]
        #n_edits.append(num_mismatch(wt, donor)) no longer have this - TODO enforce that all edits inside guide?
        edits_in_pam.append(donor[54:56]!="GG")
        n_edits_guide.append(num_mismatch(wt[33+shift:53+shift]+wt[54+shift:56+shift], donor[33:53]+donor[54:56])) 
        #used to not include edits in PAM, now it does 
        edit_pos.append(closest_edit_to_cut(wt[33+shift:56+shift], donor[33:56]))
        homopolymer.append(longest_run(donor))
    #idb["edits"] = n_edits #no longer using this
    idb["edit_in_pam"] = edits_in_pam
    idb["edits_in_guide"] = n_edits_guide 
    idb["edit_position"] = edit_pos  #TODO I'd prefer to just get 1 2 for an edit in pam - this is the case now I think
    idb["repeat_length"] = homopolymer

        
def closest_edit_to_cut(guide, edited):
    #include PAM
    assert(len(guide) == 23)
    assert(len(edited) == 23)
    edits = np.where(np.array([x for x in guide])!=np.array([y for y in edited]))[0]
    edits = edits-20
    edits = edits[edits != 0] #editing N in NGG does nothing, don't count it
    if len(edits)==0:
        return -30 #NOTE: this signifies it's not even in the guide
    else:
        return int(max(edits))
    
def longest_run(seq):
    max_c = 1
    c = 1
    for i in range(1,len(seq)):
        if seq[i] == seq[i-1]:
            c += 1
        else:
            max_c = max(max_c, c)
            c = 1
    return max_c

#CHANGE definition of 'target' to include new window.
def count_replicates(inserts_db, modify=True):
    #how many guides each target has
    replicates = []
    for i in inserts_db.index:
        replicates.append(np.sum((inserts_db.gene == inserts_db.gene[i])
                                &(inserts_db.new_window == inserts_db.new_window[i])
                                &(inserts_db.w_gene_pos == inserts_db.w_gene_pos[i]))) 
    counter = collections.Counter(replicates)
    print("Counts:",{i:counter[i] for i in range(1, max(counter)+1)})
    print("Targets 2+ replicates:", 
          sum([counter[i]/i for i in range(2, max(counter)+1)]))
    print("Targets 3+ replicates:", 
          sum([counter[i]/i for i in range(3, max(counter)+1)]))
    if modify:
        inserts_db["replicates"] = replicates
        return inserts_db
    
def within_guide_pos(insert_row, min_pos=-15, max_pos=2, ha_bp=50):
    #generally between -15 and 2 for w3
    wt_seq = np.array([ch for ch in insert_row.wt_seq])
    mod_seq = np.array([ch for ch in insert_row.donor])
    edit_pos = np.where(wt_seq != mod_seq)[0]
    edit_pos = edit_pos - ha_bp - 3 #because 0 index is 3 after cut_pos which is at 50
    if min(edit_pos) >= min_pos and max(edit_pos) <= max_pos:
        return True
    else:
        return False
    #making sure edits fall between some pos in guide

def select_guides(inserts_db, seed_size=7):
    #select guides that have < 10 bp homopolymer
    #and an edit within the seed region (-7) in the guide
    #or an edit in the PAM
    repeats = inserts_db.repeat_length < 10
    print(np.sum(~repeats), "guide have a long (10+) repeat")
    guide_edit = inserts_db.edit_position >= -seed_size
    PAM_edit = inserts_db.edit_in_pam 
    print(np.sum(PAM_edit), "guides have an edit in PAM")
    print(np.sum(~(guide_edit|PAM_edit)), "guides have edits not in PAM or seed region")
    print("Edit positions: ", sorted(collections.Counter(inserts_db.edit_position).items()))
    inserts_db["selected"] = repeats & (guide_edit | PAM_edit)
    
#what we want is to filter for unique guides
def unique_guides_filter(idb, maximize_replicates=True):
    #does this modify in place or return a new one?
    #if "target" not in idb.columns:
    #    idb["target"] = [idb.gene[i]+"_"+str(idb.w_codon_num[i]) for i in idb.index]
    #for guides that are not unique...
    #choose one to leave at random?
    #drop from target that has fewer replicates?
    #for now, just mark all but first as "not_unique"
    duplicates_list = []
    for guide in list(set(idb.guide_seq)):
        idb_guide = idb[idb.guide_seq==guide]
        if len(idb_guide) > 1:
            if maximize_replicates:
            #but this doesn't update replicate values....?
            #so sometimes guides that already have fewer replicates may be prioritized.
                best_guide = idb_guide.replicates.idxmax()
            #actually, drop ones where edit is in PAM or maximal
            elif sum(idb_guide.edit_in_pam) >= 1:
                best_guide = idb_guide.edit_in_pam.idxmax()
            else:
                best_guide = idb_guide.edit_position.idxmax()
            duplicates = list(idb_guide.index)
            duplicates.remove(best_guide)
            duplicates_list += duplicates
    idb["duplicate_guide"] = [x in duplicates_list for x in idb.index]
    return idb 

def num_targets(inserts_db):
    #assumes all windows are unique [wrong]
    if "target" in inserts_db.columns:
        return len(set(inserts_db.target))
    return len(set(inserts_db.slow))


def slow_to_fast(seq):
    if len(seq)%3 != 0:
        raise ValueError("len seq not mod 3")
    codon_seq = [seq[i*3:i*3+3] for i in range(len(seq)//3)]   
    aa_list = [codon_info.loc[c].AA for c in codon_seq]
    #print(aa_list)
    fast_list = [codon_info[codon_info.AA==aa].er.idxmax() 
                 for aa in aa_list]
    return "".join(fast_list)



def get_start_positions(seq):
    codon_seq = np.array([seq[i*3:i*3+3] for i in range(len(seq)//3)])
    start_positions = list(np.where(codon_seq == "ATG")[0])
    return start_positions #in codons

def downstream_starts(gene, codon): 
    seq = db_genes.seq[gene][codon*3+3:]
    return len(get_start_positions(seq))

def percent_loc(gene, codon):
    gene_length = (db_genes.cdsEnd[gene] - db_genes.cdsStart[gene])//3
    return codon/gene_length
