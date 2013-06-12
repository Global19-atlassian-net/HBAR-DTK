import sys
import os
from pbcore.io import FastaIO



def rev_aln_strand(aln):
    q_strand, q_s, q_e, q_l = aln[0]
    t_strand, t_s, t_e, t_l = aln[1]
    q_strand = 1 - q_strand
    q_s, q_e = q_l - q_e, q_l - q_s
    t_strand = 1 - t_strand
    t_s, t_e = t_l - t_e, t_l - t_s
    return ( (q_strand, q_s, q_e, q_l), (t_strand, t_s, t_e, t_l) )

class Overlap(object):
    """
    represent the overlap information.
    """
    def __init__(self, query_id, target_id, aln, aln_score, aln_idt, containment_tolerance):
        "the alginment data will be normalized to such the the query is always 5' -> 3'"
        self.query_id = query_id
        self.target_id = target_id
        q_strand, q_s, q_e, q_l = aln[0]
        t_strand, t_s, t_e, t_l = aln[1]
        self.aln_score = aln_score
        self.aln_idt = aln_idt
        if q_strand == 1:
            self.aln = rev_aln_strand(aln)
        else:
            self.aln = aln
        self.t_offset = 0
        self.containment_tolerance = containment_tolerance
        self.overlap_type, self.t_offset = self.get_overlap_type(self.containment_tolerance)

    def get_overlap_type(self, containment_tolerance):
        q_strand, q_s, q_e, q_l = self.aln[0]
        t_strand, t_s, t_e, t_l = self.aln[1]
        t_offset = None
        if q_s < containment_tolerance and q_e > q_l - containment_tolerance:
            """   ---------------------------->  target """
            """       --------------------->     query  """
            """   or                                    """
            """   <----------------------------  target """
            """       --------------------->     query  """
            overlap_type = "contained"      
            t_offset = - q_s
        elif t_s < containment_tolerance and t_e > t_l - containment_tolerance:
            """       --------------------->     target """
            """   ---------------------------->  query  """
            """   or                                    """
            """       <---------------------     target """
            """   ---------------------------->  query  """
            overlap_type = "contains"
            t_offset = q_s
        elif q_s <= containment_tolerance and t_l - t_e <= containment_tolerance:
            """   --------------------->         target """
            """           -------------------->  query  """
            """   or                                    """
            """   <---------------------         target """
            """           -------------------->  query  """
            overlap_type = "5p"
            t_offset = - q_s
        elif q_l - q_e <= containment_tolerance and t_s <= containment_tolerance:
            """           -------------------->  target """
            """   --------------------->         query  """
            """   or                                    """
            """           <--------------------- target """
            """   -------------------->          query  """
            overlap_type = "3p"
            t_offset = q_s
        else:
            overlap_type = "split"

        #print q_strand, q_s, q_e, q_l 
        #print t_strand, t_s, t_e, t_l
        #print overlap_type, t_offset
        return overlap_type, t_offset

class SequenceFragement(object):

    def __init__(self, id_):
        self.id_ = id_
        self.seq = None 
        self.overlaps = [] 

    def append_overlap(self, overlap):
        self.overlaps.append(overlap)

    def load_seq(self, seq):
        self.seq = seq

    def best_overlap(self, end):
        score_ovlp = []
        for ovlp in self.overlaps:
            if ovlp.overlap_type != end:
                continue
            score_ovlp.append( ovlp )
        if score_ovlp:
            score_ovlp.sort(key = lambda x: x.aln_score )
            return score_ovlp[0]
        else:
            return None 

    @property
    def best_3p_overlap(self):
        return self.best_overlap("3p")
    
    @property
    def best_5p_overlap(self):
        return self.best_overlap("5p")

def get_contained_or_split_reads_from_m4(m4_filename, permitted_error_pct, containment_tolerance):
    split_read_frgs = set() 
    contained_read_frgs = set()
    with open(m4_filename) as m4_f:
        for l in m4_f:
            l = l.strip().split()
            q_name, t_name = l[0:2]
            if q_name == t_name:
                continue
            aln_score = int(l[2])
            aln_idt = float(l[3])
            if aln_idt < 100 - permitted_error_pct:
                continue
            q_strand, q_s, q_e, q_l = ( int(x) for x in l[4:8])
            t_strand, t_s, t_e, t_l = ( int(x) for x in l[8:12])
            aln = ( (q_strand, q_s, q_e, q_l), (t_strand, t_s, t_e, t_l) )
            ovlp = Overlap(q_name, t_name, aln, aln_score, aln_idt, containment_tolerance) 
            #print " ".join(l), ovlp.overlap_type
            if ovlp.overlap_type == "split":
                if q_s < 500 or q_l - q_e < 500:
                    split_read_frgs.add( q_name )
                if t_s < 500 or t_l - t_e < 500:
                    split_read_frgs.add( t_name )
            #elif ovlp.overlap_type == "contained":
            #    contained_read_frgs.add( q_name )
            #elif ovlp.overlap_type == "contains":
            #    contained_read_frgs.add( t_name )
    with open(m4_filename) as m4_f:
        for l in m4_f:
            l = l.strip().split()
            q_name, t_name = l[0:2]
            if q_name == t_name:
                continue
            aln_score = int(l[2])
            aln_idt = float(l[3])
            if aln_idt < 100 - permitted_error_pct:
                continue
            q_strand, q_s, q_e, q_l = ( int(x) for x in l[4:8])
            t_strand, t_s, t_e, t_l = ( int(x) for x in l[8:12])
            aln = ( (q_strand, q_s, q_e, q_l), (t_strand, t_s, t_e, t_l) )
            ovlp = Overlap(q_name, t_name, aln, aln_score, aln_idt, containment_tolerance) 
            #print " ".join(l), ovlp.overlap_type
            if ovlp.overlap_type == "split":
                continue
            elif ovlp.overlap_type == "contained":
                if t_name in split_read_frgs:
                    if q_name in split_read_frgs:
                        contained_read_frgs.add( q_name )
                else:
                    contained_read_frgs.add( q_name )
            elif ovlp.overlap_type == "contains":
                if q_name in split_read_frgs:
                    if t_name in split_read_frgs:
                        contained_read_frgs.add( t_name )
                else:
                    contained_read_frgs.add( t_name )
    return split_read_frgs, contained_read_frgs

def get_unique_ovlp_reads_from_m4(m4_filename, split_reads, contained_reads, permitted_error_pct, containment_tolerance):
    read_frgs = {}
    excluded_reads = split_reads | contained_reads
    #excluded_reads =  contained_reads
    with open(m4_filename) as m4_f:
        for l in m4_f:
            l = l.strip().split()
            q_name, t_name = l[0:2]
            if q_name == t_name:
                continue
            if q_name in excluded_reads:
                continue
            if t_name in excluded_reads:
                continue

            aln_score = int(l[2])
            aln_idt = float(l[3])
            if aln_idt < 100 - permitted_error_pct:
                continue
            q_strand, q_s, q_e, q_l = ( int(x) for x in l[4:8])
            t_strand, t_s, t_e, t_l = ( int(x) for x in l[8:12])
            aln = ( (q_strand, q_s, q_e, q_l), (t_strand, t_s, t_e, t_l) )
            ovlp = Overlap(q_name, t_name, aln, aln_score, aln_idt, containment_tolerance) 
            #print " ".join(l), ovlp.overlap_type
            if ovlp.overlap_type not in ["3p", "5p"]:
                continue

            read_frgs[q_name] = read_frgs.get( q_name, SequenceFragement(q_name) )
            read_frgs[q_name].append_overlap(ovlp)

            aln1, aln2 = ovlp.aln
            aln = aln2, aln1
            read_frgs[t_name] = read_frgs.get( t_name, SequenceFragement(t_name) )
            new_ovlp = Overlap(t_name, q_name, aln, ovlp.aln_score, ovlp.aln_idt, containment_tolerance) 
            read_frgs[t_name].append_overlap(new_ovlp)

    return read_frgs

rev_map = dict(zip("ACGTacgtNn","TGCAtgcaNn"))
def rev_cmp(seq):
    return "".join([rev_map[c] for c in seq[::-1]])

def get_path(read_frgs, init_read_frg, direction="3p", visited=set()):
    reverse_orientation = { "=>":"<=", "<=":"=>" }
    #visited = set()
    cur_frg = init_read_frg
    cur_orientation = "=>"
    cum_len = 0
    path = []
    seqs = []
    c_s = 0
    c_e = len(cur_frg.seq)
    while 1:

        visited.add(cur_frg.id_)

        if direction == "3p":

            if cur_orientation == "=>":

                if cur_frg.best_3p_overlap == None:
                    seqs.append( cur_frg.seq[c_s:] )
                    path.append( (direction, cur_orientation, cur_frg.id_, c_s, c_e, 
                        next_overlap.query_id, q_strand, q_s, q_e, q_l, next_overlap.target_id, t_strand, t_s, t_e, t_l) )
                    break
                else:
                    """
                                  t_s    t_e
                                  -----------------> next_frg
                      --------------------> cur_frg
                                  q_s    q_e

                                  t_e    t_s
                                  <----------------- next_frg
                      --------------------> cur_frg
                                  q_s    q_e
                    """
                    next_overlap = cur_frg.best_3p_overlap
                    q_aln, t_aln = next_overlap.aln
                    q_strand, q_s, q_e, q_l = q_aln
                    t_strand, t_s, t_e, t_l = t_aln
                    c_e = q_s
                    seqs.append( cur_frg.seq[c_s:c_e] )
                    path.append( (direction, cur_orientation, cur_frg.id_, c_s, c_e, 
                        next_overlap.query_id, q_strand, q_s, q_e, q_l, next_overlap.target_id, t_strand, t_s, t_e, t_l) )

                    #print "1:",c_s, c_e

                    if t_strand == 0:
                        c_s = t_s
                    else:
                        c_e = t_l - t_s
                        cur_orientation = reverse_orientation[cur_orientation]

            elif cur_orientation == "<=":

                if cur_frg.best_5p_overlap == None:
                    seqs.append( rev_cmp(cur_frg.seq[:c_e]) ) 
                    path.append( (direction, cur_orientation, cur_frg.id_, c_s, c_e, 
                        next_overlap.query_id, q_strand, q_s, q_e, q_l, next_overlap.target_id, t_strand, t_s, t_e, t_l) )
                    break
                else:
                    """
                                  t_e    t_s
                                  <----------------- next_frg
                      <-------------------- cur_frg
                                  q_e    q_s

                                  t_s    t_e
                                  -----------------> next_frg
                      <-------------------- cur_frg
                                  q_e    q_s

                    """
                    next_overlap = cur_frg.best_5p_overlap
                    q_aln, t_aln = next_overlap.aln
                    q_strand, q_s, q_e, q_l = q_aln
                    t_strand, t_s, t_e, t_l = t_aln
                    c_s = q_e
                    seqs.append( rev_cmp(cur_frg.seq[c_s:c_e]) )
                    path.append( (direction, cur_orientation, cur_frg.id_, c_s, c_e, 
                        next_overlap.query_id, q_strand, q_s, q_e, q_l, next_overlap.target_id, t_strand, t_s, t_e, t_l) )

                    #print "3:",c_s, c_e

                    if t_strand == 0:
                        c_e = t_e
                    else:
                        c_s = t_l - t_e
                        cur_orientation = reverse_orientation[cur_orientation]

        if direction == "5p":

            if cur_orientation == "=>":

                if cur_frg.best_5p_overlap == None:
                    seqs.append( cur_frg.seq[:c_e] )
                    path.append( (direction, cur_orientation, cur_frg.id_, c_s, c_e, 
                        next_overlap.query_id, q_strand, q_s, q_e, q_l, next_overlap.target_id, t_strand, t_s, t_e, t_l) )
                    break
                else:
                    """
                          t_s    t_e
                    -------------> next_frg
                           --------------------> cur_frg
                          q_s    q_e

                          t_e    t_s
                    <------------- next_frg
                           --------------------> cur_frg
                          q_s    q_e
                    """
                    next_overlap = cur_frg.best_5p_overlap
                    q_aln, t_aln = next_overlap.aln
                    q_strand, q_s, q_e, q_l = q_aln
                    t_strand, t_s, t_e, t_l = t_aln
                    c_s = q_e
                    seqs.append( cur_frg.seq[c_s:c_e] )
                    path.append( (direction, cur_orientation, cur_frg.id_, c_s, c_e, 
                        next_overlap.query_id, q_strand, q_s, q_e, q_l, next_overlap.target_id, t_strand, t_s, t_e, t_l) )

                    #print "2:",c_s, c_e

                    if t_strand == 0:
                        c_e = t_e
                    else:
                        c_s = t_l - t_e
                        cur_orientation = reverse_orientation[cur_orientation]

            elif cur_orientation == "<=":

                if cur_frg.best_3p_overlap == None:
                    seqs.append( rev_cmp(cur_frg.seq[c_s:]) )
                    path.append( (direction, cur_orientation, cur_frg.id_, c_s, c_e, 
                        next_overlap.query_id, q_strand, q_s, q_e, q_l, next_overlap.target_id, t_strand, t_s, t_e, t_l) )
                    break
                else:
                    """
                               t_e    t_s
                      <----------------- next_frg
                                <-------------------- cur_frg
                               q_e    q_s

                               t_s    t_e
                      -----------------> next_frg
                                <-------------------- cur_frg
                               q_e    q_s

                    """
                    next_overlap = cur_frg.best_3p_overlap
                    q_aln, t_aln = next_overlap.aln
                    q_strand, q_s, q_e, q_l = q_aln
                    t_strand, t_s, t_e, t_l = t_aln
                    c_e = q_s
                    seqs.append( rev_cmp(cur_frg.seq[c_s:c_e]) )
                    path.append( (direction, cur_orientation, cur_frg.id_, c_s, c_e, 
                        next_overlap.query_id, q_strand, q_s, q_e, q_l, next_overlap.target_id, t_strand, t_s, t_e, t_l) )

                    #print "4:",c_s, c_e

                    if t_strand == 0:
                        c_s = t_s 
                    else:
                        c_e = t_l - t_s
                        cur_orientation = reverse_orientation[cur_orientation]



        #print direction, next_overlap.query_id, next_overlap.target_id, next_overlap.aln
        #print c_s, c_e
        #if next_overlap.aln[1][0] == 1:
        #    cur_orientation = reverse_orientation[cur_orientation]
        cur_frg = read_frgs[next_overlap.target_id]
        if cur_frg.id_ in visited:
            break

    if direction == "5p":
        seqs.reverse()

    return "".join(seqs), path

def load_seq(frgs, fasta_fn):
    f = FastaIO.FastaReader(fasta_fn)
    for r in f:
        if r.name in frgs:
            frgs[r.name].seq = r.sequence

def output_split_reads(split_reads, contained_reads, fasta_fn, out_fn):
    f = FastaIO.FastaReader(fasta_fn)
    included = split_reads - contained_reads
    with open(out_fn, "w") as out:
        for r in f:
            if r.name in included:
                print >> out, ">"+r.name
                print >> out, r.sequence


if __name__ == "__main__":
    #blasr preads.fa preads.fa -maxScore -150 -m 4 -nproc 32 -bestn 64 -nCandidates 128 -maxLCPLength 25 -minMatch 24 -scoreMatrix "-1 10 10 10 10  10 -1 10 10 10  10 10 -1 10 10  10 10 10 -1 10  10 10 10 10 10" -indel 10 -noSplitSubreads -out pr_pr_strigent.m4
    import sys

    pread_fn = sys.argv[1]
    aln_m4_fn = sys.argv[2]
    output_prefix = sys.argv[3]

    split_reads, contained_reads =  get_contained_or_split_reads_from_m4(aln_m4_fn, 4, 100)
    uniq_ovlp_reads = get_unique_ovlp_reads_from_m4( aln_m4_fn, split_reads, contained_reads, 4, 100)
    
    load_seq(uniq_ovlp_reads, pread_fn)
    output_split_reads(split_reads, contained_reads, pread_fn, output_prefix+"_split.fa")

    seeds = []
    all_reads = set()
    for read_id, read_frg in uniq_ovlp_reads.items():
        all_reads.add(read_id)
        #print ctg_id, ctg_frg.best_5p_overlap, ctg_frg.best_3p_overlap
        if read_frg.best_5p_overlap == None and read_frg.best_3p_overlap != None:
            seeds.append( (len(read_frg.seq), read_id, read_frg, "3p") )
        if read_frg.best_3p_overlap == None and read_frg.best_5p_overlap != None:
            seeds.append( (len(read_frg.seq), read_id, read_frg, "5p") )

    seeds.sort(reverse=True)
    visited = set()
    ctg_id = 0
    ctg_layout_fn = output_prefix + ".lay"
    ctg_layout = open(ctg_layout_fn, "w")
    with open(output_prefix + ".fa","w") as f:
        for read_len, read_id, read_frg, direction in seeds:
            if read_id not in visited:
                seq, path = get_path(uniq_ovlp_reads, read_frg, direction=direction, visited = visited)
                if len(seq) < 200:
                    continue
                print >>f, ">ctg_%06d" % ctg_id
                print >>f, seq
                print >> ctg_layout, ">ctg_%06d" % ctg_id
                for p in path:
                    print >> ctg_layout, " ".join([str(c) for c in p])
                ctg_id += 1

        for read_id in all_reads - visited:
            if read_id not in visited:
                read_frg = uniq_ovlp_reads[read_id]
                seq, path = get_path(uniq_ovlp_reads, read_frg, direction="3p", visited = visited)
                if len(seq) < 200:
                    continue
                print >>f, ">ctg_%06d_s" % ctg_id
                print >>f, seq
                print >> ctg_layout, ">ctg_%06d_s" % ctg_id
                for p in path:
                    print >> ctg_layout, " ".join([str(c) for c in p])
                ctg_id += 1
        print >> ctg_layout, "#", all_reads - visited, len(uniq_ovlp_reads)
    ctg_layout.close()

