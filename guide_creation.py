import pandas as pd
import numpy as np
import re
from ast import literal_eval
import collections
import os

#load_files
db_genes = pd.read_csv("data/db_genes.csv", sep=",", index_col=0, header=0)
db_genes = db_genes.where(pd.notnull(db_genes), None) 
codon_info = pd.read_pickle("data/codon_info.pkl")

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

Outputs:
- inserts_db: Pandas dataframe containing the donor-guide CRISPEY insert sequences along with relevant information 
'''
def make_library(genes, positions, original_windows, deletion=True, w=3, do_not_filter=False):
    assert(len(genes)==len(positions))
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
        if gene_pos_to_chrom_pos(gene,pos,w=w, throw_error=False) is None:
            #skip if window not fully in one exon,
            #make inserts would throw error.
            continue
        idb = make_inserts(gene, pos, original_windows[i], w=w, deletion=deletion)
        #sanity_checks(idb, w=w) for speed, will do later
        # add stdev info, too.
        idb_list.append(idb)
    inserts_db = pd.concat(idb_list)
    inserts_db = inserts_db.reset_index()
    inserts_db = inserts_db.drop("index", axis=1) 
    inserts_db["target"] = [inserts_db.gene[i]+"_"+str(inserts_db.w_gene_pos[i])
                          for i in inserts_db.index]
    if do_not_filter:
        print("DID NOT FILTER")
        return inserts_db
    #instead of dropping mismatches/etc., lets keep them?
    print("")
    print("Before filtering for off-target matches:", len(inserts_db), 
          "guides for", len(set(list(inserts_db.target))), "targets")
    inserts_db = filter_guides(inserts_db)
    inserts_db = inserts_db[inserts_db.n_off_targets == 0]
    print("running sanity checks")
    sanity_checks(inserts_db, deletion=deletion)
    print("")
    print("After filtering for off-target matches:", len(inserts_db), 
          "guides for", len(set(list(inserts_db.target))), "targets")
    add_edits_info(inserts_db)
    inserts_db = count_replicates(inserts_db)
    select_guides(inserts_db)
    #TODO: do I want to ensure that guides are unique? nah
    inserts_db = inserts_db[inserts_db.selected]
    print("")
    print("After selecting for repeats/edits in guide:", len(inserts_db), 
          "guides for", len(set(list(inserts_db.target))), "targets")
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
def make_inserts(gene, gene_pos, original_window, w=3, h_arm_length=50, alt_fast=True, deletion=True):
    #gene_pos, w, and h_arm_length are in nucleotides
    #if deletion is False, change a window of w/3 codons to the fastest possible codons
    #if deletion is True, delete a window of w bp
    if deletion == False:
        if gene_pos % 3 != 0:
            raise ValueError('Gene pos is not in frame')
        if w % 3 != 0:
            raise ValueError('Window is not mod 3 - out of frame')
    assert(len(original_window) == w)
    constant_5ph = "GAGTTACTGTCTGTTTTCCT"
    constant_rt = "AGGAAACCCGTTTCTTCTGACGTAAGGGTGCGCA"
    constant_3ph = "GTTTCAGAGCTATGCTGGAA"
    chrom, chrom_pos, gene_strand = gene_pos_to_chrom_pos(gene,gene_pos,w=w, w_in_bp=True) 
    guide_list = get_all_guides(chrom, chrom_pos, w=int(np.ceil(w//3))) #w for get all guides is in codons....?
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
    new_window = "" #if deletion
    if not deletion: #change to a new window
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
    inserts_db["new_window"] = new_window
    wt_donor_list = []
    mod_donor_list = []
    insert_list = []
    for i in inserts_db.index:
        #donor direction --> should match guide direction ??
        cp = inserts_db.cut_pos[i]
        guide_strand = inserts_db.guide_strand[i]
        guide_seq = inserts_db.guide_seq[i]
        if not deletion:
            #if no indels, just substitution
            donor = get_seq(chrom, cp-h_arm_length, cp+h_arm_length)
            loc_w = h_arm_length + chrom_pos - cp
        if deletion:
            #depending on if chrom_pos is before/after cut_pos
            shift = max(min(cp - chrom_pos, w), 0) #so between del_size (w) and 0
            donor = get_seq(chrom, cp-h_arm_length-shift, cp+h_arm_length+w-shift) #+ strand
            loc_w = h_arm_length + chrom_pos - cp + shift    
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
    return inserts_db

'''Returns all possible guide sequences that are max_dist away from a starting 
position, checking both strands'''
def get_all_guides(chrom,start,w=3,max_dist=30): 
    #w and max_dist are in bp
    end = start+w
    guide_length = 20
    back_dist = max_dist-w-4 #distance to search back 
    #(4-> 'cut' 1 2 3 N | G G --> to limit where to look for GG)
    fw_dist = max_dist-w+6 #distance to search fw
    #(4-> 'cut' 1 2 3 N G G | --> to limit where to look for GG)
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
def sanity_checks(idb, additional_checks=False, deletion=False): #assumes 50 bp arms, cut pos in middle
    for i in idb.index:
        w = idb.w_size[i] #IN BP
        assert(len(idb.donor[i])==100)
        shift = 0
        del_size = 0
        guide_s = 1
        if idb.guide_strand[i]=="-":
            guide_s = -1
        if deletion:
            shift = min(max(idb.cut_pos[i]-idb.w_chrom_pos[i], 0), w)
            del_size = w
            if idb.guide_strand[i]=="-":
                shift = del_size-shift #I think?
        assert(len(idb.wt_seq[i])+del_size==100)
        assert(idb.wt_seq[i][33+shift:53+shift] == idb.guide_seq[i]), str(idb.gene[i]) + "_" + str(idb.w_codon_num[i])
        assert(idb.wt_seq[i][54+shift:56+shift]=="GG") #PAM assertion 
        assert(idb.original_window[i]==idb.old_window[i]), i
        assert(len(idb.insert_seq[i])==194)
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
            assert(donor_on_strand[50+idb.w_chrom_pos[i]-idb.cut_pos[i]:50+idb.w_chrom_pos[i]-idb.cut_pos[i]+w*3]==w_on_strand)
            nm_window = num_mismatch(idb.old_window[i], idb.new_window[i])
            nm_donor = num_mismatch(idb.original_window[i], idb.donor[i])  
            assert(nm_window==nm_donor)
            assert(len(idb.new_window[i])==w)
        else:
            insert_pos = 50 + (idb.w_chrom_pos[i] - idb.cut_pos[i])*guide_s + shift
            if guide_s == -1:
                insert_pos = insert_pos - del_size
            insert_seq = idb.old_window[i]
            assert(len(insert_seq)==del_size)
            if idb.guide_strand[i] != idb.gene_strand[i]:
                insert_seq = rev_complement(insert_seq)
            assert_msg = idb.gene[i] + idb.gene_strand[i] + idb.guide_strand[i] + " "
            assert_msg += insert_seq + " " + str(shift) +  " " + str(insert_pos) + " "+ idb.donor[i] + " " + idb.wt_seq[i]
            assert(idb.donor[i][:insert_pos]+insert_seq+idb.donor[i][insert_pos:]==idb.wt_seq[i]), assert_msg
            assert(idb.wt_seq[i][insert_pos:insert_pos+del_size]==insert_seq), assert_msg
        
        if additional_checks:
            msg = [idb.gene[i], idb.w_codon_num[i], idb.w_size[i]]
            #check that window does not cross introns/ends
            assert(not intersecting_intron(idb.gene[i], idb.w_chrom_pos[i], idb.w_chrom_pos[i]+w)), msg
            #check that slow is indeed correct
            assert(get_seq_in_gene(idb.gene[i], idb.w_gene_pos[i], idb.w_gene_pos[i]+w) == idb.original_window[i])
            if idb.gene_strand[i] == '-':
                assert(rev_complement(get_seq(idb.chrom[i], idb.w_chrom_pos[i], idb.w_chrom_pos[i]+w)) == idb.original_window[i])
            else:
                assert(get_seq(idb.chrom[i], idb.w_chrom_pos[i], idb.w_chrom_pos[i]+w) == idb.original_window[i])  
            #a change was made!
            assert(idb.old_window[i] != idb.new_window[i])
            assert(idb.wt_seq[i] != idb.donor[i])
            #check that amino acid sequence slow, fast is the same
            if not deletion:
                assert(amino_acid_seq(idb.old_window[i]) == amino_acid_seq(idb.new_window[i]))
                #fast is faster than slow?  #ah but for neutral [1 codon change] it's not actually faster .... oops
                if idb.w_size[i] != 3: #the following is not true for single codon changes (where I used alt_fast)
                    assert(slow_to_fast(idb.old_window[i]) == idb.new_window[i]) #
                assert(idb.w_codon_num[i] * 3 == idb.w_gene_pos[i])
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
def get_seq(chrom,start,end):
    with open("yeast_genes/" + chrom + ".fa", "r") as f:
        lines = f.readlines()
        seq = "".join(lines[1:])
        seq = "".join(seq.split("\n"))
    return seq[start:end]


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
def make_fastq_file(fname, guides_dict):
    #print(guides_dict)
    with open(fname, 'w') as f:
        for key in guides_dict:
            f.write(">"+key+"\n")
            f.write(guides_dict[key]+"\n")
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
    exon_starts = literal_eval(db_genes.exonStarts[gene])
    exon_ends = literal_eval(db_genes.exonEnds[gene])
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
    exon_starts = literal_eval(db_genes.exonStarts[gene])
    exon_ends = literal_eval(db_genes.exonEnds[gene])
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
    exon_starts = literal_eval(db_genes.exonStarts[gene])
    exon_ends = literal_eval(db_genes.exonEnds[gene])
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
    exon_starts = np.array(literal_eval(db_genes.exonStarts[gene]))
    exon_ends = np.array(literal_eval(db_genes.exonEnds[gene]))
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
    
def filter_guides(inserts_db, guide_col="guide_seq"):
    #filter guides for those that have no off-targets (<= 3 mismatches) 
    #add 'n_off_targets' to 'filter_guides'
    save_f = "tmp/inserts.fa"
    #TODO
    make_fastq_file(save_f, {str(i):inserts_db[guide_col][i] for i in inserts_db.index})
    #%cd /mnt/lareaulab/shelen/slow_codons_2021/data/
    os.system("bowtie -f -v 3 -a yeast_genes/sacCer_index tmp/inserts.fa > tmp/bowtie_out.csv")
    #%cd /mnt/lareaulab/shelen/slow_codons_2021/code/
    load_f = "tmp/bowtie_out.csv"
    bowtie_info = pd.read_csv(load_f, sep="\t", header=None, names=["guide_id","strand","chr","chr_pos","seq","aln","?","mismatch"])
    #maybe there's a lot of guides that only align to one spot, but have multiple targets
    inserts_db["n_off_targets"] = [sum(bowtie_info.guide_id == i)-1 for i in inserts_db.index]
    return inserts_db
    #guides_to_use = [seq for seq in inserts_db[guide_col] if sum(bowtie_info.seq==seq) == 1]
    #return guides_to_use #sequences of guides

def add_edits_info(idb): #assumes 50 bp arms
    #add info about edits to db
    #how many edits away from wt
    #how many edits in guide
    #position of edit closest to cut site?
    # -17 .... -1 |cut| 0 1 2   (indexing scheme)
    #change indexing scheme to -
    #longest homopolymer in [mutated] donor 
    n_edits = []
    n_edits_guide = []
    edits_in_pam = []
    edit_pos = []
    homopolymer = []
    for i in idb.index:
        wt = idb.wt_seq[i]
        donor = idb.donor[i]
        n_edits.append(num_mismatch(wt, donor))
        edits_in_pam.append(donor[54:56]!="GG")
        n_edits_guide.append(num_mismatch(wt[33:53], donor[33:53])) #does not include edits in PAM
        edit_pos.append(closest_edit_to_cut(wt[33:56], donor[33:56]))
        homopolymer.append(longest_run(donor))
    idb["edits"] = n_edits
    idb["edit_in_pam"] = edits_in_pam
    idb["edits_in_guide"] = n_edits_guide 
    idb["edit_position"] = edit_pos  #TODO I'd prefer to just get 1 2 for an edit in pam
    idb["repeat_length"] = homopolymer
    
        
        
def closest_edit_to_cut(guide, edited):
    #include PAM
    assert(len(guide) == 23)
    assert(len(edited) == 23)
    edits = np.where(np.array([x for x in guide])!=np.array([y for y in edited]))[0]
    edits = edits-20
    edits = edits[edits != 0]
    if len(edits)==0:
        return None
    else:
        return max(edits)
    
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

def count_replicates(inserts_db, modify=True):
    #how many guides each target has
    replicates = []
    for i in inserts_db.index:
        replicates.append(np.sum((inserts_db.gene == inserts_db.gene[i])&(inserts_db.w_gene_pos == inserts_db.w_gene_pos[i]))) 
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

def select_guides(inserts_db):
    #select guides that have < 10 bp homopolymer
    #and an edit within the seed region (-7) in the guide
    #or an edit in the PAM
    repeats = inserts_db.repeat_length < 10
    guide_edit = inserts_db.edit_position >= -7 
    PAM_edit = inserts_db.edit_in_pam
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

