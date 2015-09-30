#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: tgibbons (tgibbons@umd.edu)
"""

# Standard Python libraries
from sys import stderr
import sys
import argparse

# Third-party libraries
import networkx as nx

def main(argv=None):
    """Where the magic happens!

    The main() function coordinates calls to all of the other functions in this
    program in the hope that, by their powers combined, useful work will be
    done.

    Args:
        None

    Returns:
        An exit status (hopefully 0)
    """
    if argv is None:
        argv = sys.argv

    args = get_parsed_args()

    # Check for sane runtime arguments and data
    len_cols = argument_error_handling(args.blast)
    if not len_cols:
        min_cov = False
    else:
        min_cov = args.min_cov

    graph = nx.Graph()

    if args.combine:
        add_greedily_combined_edges(graph=graph,
                                    blast_handle=args.blast,
                                    min_cov=min_cov)
    else:
        add_top_hit_edges(graph=graph,
                          blast_handle=args.blast,
                          min_cov=min_cov)

    if args.norm:
        avgs = compute_organism_averages(graph=graph, idchar=args.idchar)
        normalize_graph(graph=graph, avgs=avgs, idchar=args.idchar)

    print_abc_file(graph=graph, abc=args.abc)

def add_greedily_combined_edges(graph, blast_handle, min_cov=False):
    """Construct a graph by combining non-overlapping hits between each pair of sequences"""
    crnt_qry = None  # Current query sequence ID
    crnt_ref = None  # Current reference sequence ID
    hit_stack = []  # Stack of tuples with bit scores and coordinates

    for line in blast_handle:
        temp = line.strip().split()
        if not temp:  # skip blank lines
            continue
        elif temp[0][0] == "#":  # skip comment lines
            continue
        elif temp[0] == temp[1]:  # skip self-hits
            continue

        qry_id = str(temp[0])
        ref_id = str(temp[1])

        # This stack of commands just generates a tuple of alignment info
        bit = float(temp[11])
        id_frac = float(temp[2])/100
        aln_len = int(temp[3])
        id_sum = id_frac*aln_len
        qry_coords = sorted([int(temp[6]), int(temp[7])])
        ref_coords = sorted([int(temp[8]), int(temp[9])])
        aln_info = tuple([bit, id_sum] + qry_coords + ref_coords)

        # Keep adding alignment info to stack if the seqs remain the same
        if qry_id == crnt_qry and ref_id == crnt_ref:
            hit_stack.append(aln_info)

        # Catch first pair of seqs and initialize per-pair variables
        elif not (crnt_qry and crnt_ref):
            crnt_qry = qry_id
            crnt_ref = ref_id
            hit_stack = [aln_info]
            if min_cov:
                min_len = min(int(temp[12]), int(temp[13]))
            else:
                min_len = False

        # Process stack of hits if a seq changes
        else:
            greedily_combine_hits_into_edge(
                graph, crnt_qry, crnt_ref, hit_stack, min_len, min_cov)

            # Reset per-pair variables for next seq pair
            crnt_qry = qry_id
            crnt_ref = ref_id
            hit_stack = [aln_info]
            if min_cov:
                min_len = min(int(temp[12]), int(temp[13]))
            else:
                min_len = False

    # Process final seq pair
    if len(hit_stack) > 0:
        greedily_combine_hits_into_edge(
            graph, crnt_qry, crnt_ref, hit_stack, min_len, min_cov)

def greedily_combine_hits_into_edge(
        graph, qry_id, ref_id, hit_stack, min_len, min_cov=False):
    """Greedily combine non-overlapping hits and add resulting edge to graph"""
    # Remove information from first alignment and use it to initialize variables
    bit, id_sum, qbeg1, qend1, rbeg1, rend1 = hit_stack.pop(0)
    nlp_bit    = bit
    nlp_id_sum = id_sum
    nlp_coords = [[qbeg1, qend1, rbeg1, rend1]]

    # If additional alignments exist, check if they overlap and add them
    try:
        for aln_info in hit_stack:
            bit, id_sum, qbeg1, qend1, rbeg1, rend1 = aln_info  # Candidate hit
            for nlp_coord in nlp_coords:
                qbeg2, qend2, rbeg2, rend2 = nlp_coord  # Coordinates for accepted hit
                if qbeg1 > qend2 or qend1 < qbeg2:  # No conflicting overlaps in query
                    if rbeg1 > rend2 or rend1 < rbeg2:  # No conflicting overlaps in ref
                        nlp_bit    += bit     # Add bit score to seq pair
                        nlp_id_sum += id_sum  # Add identical character count to seq pair

    # Move along gracefully if there's only a single alignment for the seq pair
    except IndexError:
        pass

    # Make sure coverage (% identical bases) exceeds user-defined threshold
    if min_cov:
        nlp_cov = nlp_id_sum / min_len
        if nlp_cov < min_cov:
            return

    # G[A][B] != G[B][A] in nx Graphs, so add edges in node-sorted order
    id1, id2 = sorted([qry_id, ref_id])
    try:
        if nlp_bit > graph[id1][id2]['bit']:
            graph[id1][id2]['bit'] = nlp_bit
    except KeyError:
        graph.add_edge(id1, id2, bit=nlp_bit)

def add_top_hit_edges(graph, blast_handle, min_cov=False):
    """Construct a graph from the top BLAST hits"""
    for line in blast_handle:
        temp = line.strip().split()
        if not temp:  # skip blank lines
            continue
        elif temp[0][0] == "#":  # skip comment lines
            continue
        elif temp[0] == temp[1]:  # skip self-hits
            continue

        qry_id = str(temp[0])
        ref_id = str(temp[1])

        bit = float(temp[11])
        id_frac = float(temp[2])/100
        aln_len = int(temp[3])
        id_sum = id_frac*aln_len

        if min_cov:
            min_len = min(int(temp[12]), int(temp[13]))
            cov = id_sum/min_len
            if cov < min_cov:
                continue

        id1, id2 = sorted([qry_id, ref_id])
        try:
            if bit > graph[id1][id2]['bit']:
                graph[id1][id2]['bit'] = bit
        except KeyError:
            graph.add_edge(id1, id2, bit=bit)

def compute_organism_averages(graph, idchar='|'):
    """Compute symmetrical average scores between and within each pair of organisms

    :param graph: NetworkX Graph() containing BLAST graph weighted with bit scores
    :param idchar: Character used to separate organism identifier from the rest
        of the sequence identifier
    :return avgs: Dictionary containing the total number of edges between each
        pair of organisms, the cumulative sum of each metric, and the average
        score for each metric (one node per organism, one edge per pair)
    """
    avgs = dict()

    for id1, id2, edata in graph.edges(data=True):
        org1 = id1.split(idchar)[0]
        org2 = id2.split(idchar)[0]
        orgA, orgB = sorted([org1, org2])

        try:
            avgs[orgA][orgB]['cnt'] += 1
            avgs[orgA][orgB]['sum'] += edata['bit']
        except KeyError:
            try:
                avgs[orgA][orgB] = {'cnt': 1, 'sum': edata['bit']}
            except KeyError:
                avgs[orgA] = {orgB: {'cnt': 1, 'sum': edata['bit']}}

    for orgA, orgA_data in avgs.iteritems():
        for orgB, inter_org_stats in orgA_data.iteritems():
            avgs[orgA][orgB]['avg'] = inter_org_stats['sum']/inter_org_stats['cnt']

    return avgs

def normalize_graph(graph, avgs, idchar):
    """Normalize graph using inter- & intra-organism averages

    This function inflates edge weights between less similar organisms, and
    deflates edge weights between more similar organism, converting all edge
    weights into dimensionless multiples of the overall average graph edge
    weight.

    :param graph: NetworkX Graph() containing BLAST graph weighted with bit scores
    :param avgs: Dictionary containing the total number of edges between each
        pair of organisms, the cumulative sum of each metric, and the average
        score for each metric (one node per organism, one edge per pair)
    """
    for id1, id2, edata in graph.edges(data=True):
        org1 = id1.split(idchar)[0]
        org2 = id2.split(idchar)[0]
        orgA, orgB = sorted([org1, org2])
        graph[id1][id2]['bit'] /= avgs[orgA][orgB]['avg']

def print_abc_file(graph, abc):
    """Print graph file in "abc" format

    The "abc" format is very simple. The first two columns are a pair of node
    IDs, and the third column is the weight of an edge between those nodes.
    The format implies a direction for each edge and supports both unweighted
    edges and multiple edges between a pair of nodes, but node of those
    features are used here.

    :param graph: Dictionary containing BLAST graph weighted with bit scores
    :param abc: Open output file handle for "abc" graph
    """
    for id1, id2, edata in sorted(graph.edges(data=True)):
        out_line = '\t'.join([id1, id2, str(edata['bit'])])+'\n'
        abc.write(out_line)

def print_connected_component_fasta_files(graph, fasta_handle, out_pref):
    """Split FASTA file by connected component

    TODO: Make this function accessible from the main interface because it can
    be very useful for downstream analyses. The FASTA file could also be mined
    for sequence lengths in cases where they are not included in the BLAST file.
    """
    from Bio import SeqIO
    fasta = SeqIO.to_dict(SeqIO.parse(fasta_handle, 'fasta'))
    w = len(str(len(nx.connected_components(graph))))
    cmp_cnt = 0
    for comp in nx.connected_components(graph):
        cmp_hdl = open(out_pref+"_comp"+str(cmp_cnt).zfill(w)+".fasta", 'w')
        for sid in comp:
            try:
                cmp_hdl.write(">{0}\n{1}\n".format(sid, fasta[sid].seq))
            except KeyError:
                stderr.write("{0} not found in FASTA file\n".format(sid))
        cmp_hdl.close()
        cmp_cnt += 1

def get_parsed_args():
    """Parse the command line arguments

    Parses command line arguments using the argparse package, which is a
    standard Python module starting with version 2.7.

    Args:
        None, argparse fetches them from user input

    Returns:
        args: An argparse.Namespace object containing the parsed arguments
    """
    parser = argparse.ArgumentParser(
        description=("Generate a set of graphs from a tab-delimited BLASTP or "
        "BLASTN file that contains the query and subject IDs in the first two "
        "columns and the bit score in the twelfth column, unless otherwise "
        "specified"))

    parser.add_argument('--idchar',
                        dest='idchar',
                        action='store',
                        default='|',
                        help=("The character used to separate the organism ID"
                        "from the rest of the sequence header [def='|']"))

    parser.add_argument('--norm',
                        dest='norm',
                        action='store_true',
                        default=False,
                        help=("Normalize edge weights using inter-organism "
                        "averages [def=False]"))

    parser.add_argument('--min_cov',
                        dest='min_cov',
                        type=float,
                        action='store',
                        default=False,
                        help=("Minimum fraction (e.g. 0.5) of characters in "
                        "the shorter sequence that must be aligned to "
                        "identical characters in the longer sequence "
                        "[def=False]"))

    parser.add_argument('--combine',
                        dest='combine',
                        action='store_true',
                        default=False,
                        help=("Greedily combine bit scores from "
                        "non-overlapping hits [def=False]"))

    parser.add_argument('blast',
                        nargs='?',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        help=("Tab-delimited BLAST file containing additional "
                        "columns for the query and sequence lengths "
                        "(comment lines are okay) [def=stdin]"))

    parser.add_argument('abc',
                        nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="'abc' output graph file [def=stdout]")

    args = parser.parse_args()

    return args

def argument_error_handling(blast):
    """Check to make sure user has specified sane arguments

    Check BLAST file for standard bit score column and non-standard sequence
    length columns
    """
    len_cols = False

    for line in blast:
        if line[0] == '#':  # Skip comment lines
            if line[2:8] == 'Fields':
                header = line.strip().split(', ')

                # Verify standard columns using BLAST header
                standard_warning = ("Warning: It appears from the headers that "
                    "the {0} is not listed in the {1} column. Execution will "
                    "continue, although results may be incorrect.\n")

                # Query sequence column
                try:
                    if not [i for i, word in enumerate(header)
                            if word.endswith('query id')][0] == 0:
                        sys.stderr.write(standard_warning.format(
                            'query sequence ID', 'first'))
                except IndexError:
                    sys.stderr.write(standard_warning.format(
                        'query sequence ID', 'first'))

                # Subject sequence column
                try:
                    if not header.index('subject id') == 1:
                        sys.stderr.write(standard_warning.format(
                            'subject sequence ID', 'second'))
                except ValueError:
                    sys.stderr.write(standard_warning.format(
                        'subject sequence ID', 'second'))

                # Percent identity column
                try:
                    if not header.index('% identity') == 2:
                        sys.stderr.write(standard_warning.format(
                            '% identity', 'third'))
                except ValueError:
                    sys.stderr.write(standard_warning.format(
                        '% identity', 'third'))

                # Alignment length column
                try:
                    if not header.index('alignment length') == 3:
                        sys.stderr.write(standard_warning.format(
                            'alignment length', 'fourth'))
                except ValueError:
                    sys.stderr.write(standard_warning.format(
                        'alignment length', 'fourth'))

                # Column with index for starting position in query sequence
                try:
                    if not header.index('q. start') == 6:
                        sys.stderr.write(standard_warning.format(
                            'index for the start site in the query sequence',
                            'seventh'))
                except ValueError:
                    sys.stderr.write(standard_warning.format(
                        'index for the start site in the query sequence',
                        'seventh'))

                # Column with index for ending position in query sequence
                try:
                    if not header.index('q. end') == 7:
                        sys.stderr.write(standard_warning.format(
                            'index for the end site in the query sequence',
                            'eighth'))
                except ValueError:
                    sys.stderr.write(standard_warning.format(
                        'index for the end site in the query sequence',
                        'eighth'))

                # Column with index for starting position in subject sequence
                try:
                    if not header.index('s. start') == 8:
                        sys.stderr.write(standard_warning.format(
                            'index for the start site in the subject sequence',
                            'ninth'))
                except ValueError:
                    sys.stderr.write(standard_warning.format(
                        'index for the start site in the subject sequence',
                        'ninth'))

                # Column with index for ending position in subject sequence
                try:
                    if not header.index('s. end') == 9:
                        sys.stderr.write(standard_warning.format(
                            'index for the end site in the subject sequence',
                            'tenth'))
                except ValueError:
                    sys.stderr.write(standard_warning.format(
                        'index for the end site in the subject sequence',
                        'tenth'))

                # Bit score column
                try:
                    if not header.index('bit score') == 11:
                        sys.stderr.write(standard_warning.format(
                            'bit score', 'twelfth'))
                except ValueError:
                    sys.stderr.write(standard_warning.format(
                        'bit score', 'twelfth'))

                # Check for non-standard columns containing sequence lengths
                len_col_warning =  ("Warning: It appears from the headers that "
                    "the BLAST file does not contain columns for the query and "
                    "subject sequence lengths. These are necessary for the "
                    "percent match filter that requires a minimum fraction of "
                    "the characters in the shorter sequence to be identically "
                    "matched within a(n) (set of) alignment(s). The filter "
                    "will therefore be disabled. To re-enable this feature, "
                    "rerun BLAST with option -outfmt '<6/7> std qlen slen'.\n")

                try:
                    qlen_col = header.index('query length')
                    slen_col = header.index('subject length')
                    if qlen_col == 12 and slen_col == 13:
                        len_cols = True
                    else:
                        sys.stderr.write(standard_warning)
                except ValueError:
                    sys.stderr.write(len_col_warning)

                # Other tests are unnecessary when headers are present
                break

            else:
                continue

        elif not line.strip():  # Skip blank lines
            continue

        temp = line.strip().split()

        if len(temp) < 12 and bscol == 12:
            sys.stderr.write("Error: BLAST file does not appear to contain "
            "all twelve standard tab-delimited columns. It is difficult to be "
            "more specific without having header lines.\n")
            raise
        elif len(temp) < 14:
            sys.stderr.write("Warning: BLAST file does not appear to contain "
            "all sequence length information. In order to filter alignments by "
            "coverage, BLAST must be re-executed using option -outfmt '<6/7> "
            "std qlen slen'. Execution will continue without coverage "
            "filter...\n")
        else:
            len_cols = True

        break

    blast.seek(0)

    return len_cols

if __name__ == "__main__":
    sys.exit(main())
