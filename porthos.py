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

    if args.allowed:
        graph = initialize_graph(allowed=args.allowed)
        allowed = True
    else:
        graph = nx.Graph()
        allowed = False

    if args.disallowed:
        disallowed = set(line.strip() for line in open(args.disallowed))
        graph.remove_nodes_from(disallowed)
    else:
        graph = nx.Graph()
        disallowed = False

    if args.combine:
        add_greedily_combined_edges(graph=graph,
                                    blast_handle=args.blast,
                                    bscol=args.bscol-1,
                                    min_cov=args.min_cov,
                                    allowed=allowed,
                                    disallowed=disallowed)
    else:
        add_top_hit_edges(graph=graph,
                          blast_handle=args.blast,
                          bscol=args.bscol-1,
                          min_cov=args.min_cov,
                          allowed=allowed,
                          disallowed=disallowed)

    if args.norm:
        avgs = compute_organism_averages(graph=graph, idchar=args.idchar)
        normalize_graph(graph=graph, avgs=avgs, idchar=args.idchar)

    print_abc_file(graph=graph, abc=args.abc)


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
        description="Generate a set of graphs from a tab-delimited BLASTP " +
                    "or BLASTN file that contains the query and subject IDs " +
                    "in the first two columns and the bit score in the " +
                    "twelfth column, unless otherwise specified")

    parser.add_argument('--idchar',
                        dest='idchar',
                        action='store',
                        default='|',
                        help="The character used to separate the organism " +
                             "ID from the rest of the sequence header " +
                             "[def='|']")

    parser.add_argument('--bscol',
                        dest='bscol',
                        action='store',
                        default=12,
                        help="One-indexed column containing the bit scores " +
                             "[def=12]")

    parser.add_argument('--norm',
                        dest='norm',
                        action='store_true',
                        default=False,
                        help="Normalize edge weights using inter-organism " +
                             "averages [def=False]")

    parser.add_argument('--min_cov',
                        dest='min_cov',
                        type=float,
                        action='store',
                        default=0.5,
                        help="Minimum fraction of characters in the shorter " +
                             "sequence that must be aligned to identical " +
                             "characters in the longer sequence [def=0.5]")

    parser.add_argument('--combine',
                        dest='combine',
                        action='store_true',
                        default=False,
                        help="Greedily combine bit scores from " +
                             "non-overlapping hits [def=False]")

    parser.add_argument('--allowed',
                        dest='allowed',
                        action='store',
                        default=False,
                        help="File containing a list of acceptable nodes")

    parser.add_argument('--disallowed',
                        dest='disallowed',
                        action='store',
                        default=False,
                        help="File containing a list of forbidden nodes")

    parser.add_argument('blast',
                        nargs='?',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        help="Tab-delimited BLAST file containing " +
                             "additional columns for the query and sequence " +
                             " lengths (comment lines are okay) [def=stdin]")

    parser.add_argument('abc',
                        nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="'abc' output graph file [def=stdout]")

    args = parser.parse_args()

    return args


def initialize_graph(allowed=False):
    """
    """
    graph = nx.Graph()

    if allowed:
        nodes = set()
        for line in open(allowed):
            node = str(line.strip())
            if not node:  # skip blank lines
                continue
            nodes.add(node)

        graph.add_nodes_from(list(nodes))
        del nodes

    return graph


def add_greedily_combined_edges(graph, blast_handle, bscol=11, min_cov=0.5,
              allowed=False, disallowed=False):
    """Construct a graph from the top BLAST hits
    """
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

        if allowed:
            if not (graph.has_node(qry_id) and graph.has_node(ref_id)):
                continue

        elif disallowed:
            if qry_id in disallowed or ref_id in disallowed:
                continue

        # This stack of commands just generates a tuple of alignment info
        bit = float(temp[bscol])
        id_frac = float(temp[2])/100
        aln_len = int(temp[3])
        id_sum = id_frac*aln_len
        qry_coords = sorted([int(temp[6]), int(temp[7])])
        ref_coords = sorted([int(temp[8]), int(temp[9])])
        aln_info = tuple([bit, id_sum] + qry_coords + ref_coords)

        # Keep adding alignment info to stack if the seqs are the same
        if qry_id == crnt_qry and ref_id == crnt_ref:
            hit_stack.append(aln_info)
            continue

        # If seqs change, process hit
        else:
            # ...unless it's the first hit in the filed
            if len(hit_stack) > 0:
                add_greedily_combined_hits(
                    graph, crnt_qry, crnt_ref, hit_stack, min_len, min_cov)

            # Set variables used for upcoming (set of) hit(s)
            crnt_qry = qry_id
            crnt_ref = ref_id
            min_len = min(int(temp[12]), int(temp[13]))
            hit_stack = [aln_info]

    # Print final edge
    if len(hit_stack) > 0:
        add_greedily_combined_hits(
            graph, crnt_qry, crnt_ref, hit_stack, min_len, min_cov)


def add_greedily_combined_hits(
        graph, qry_id, ref_id, hit_stack, min_len, min_cov):
    """Greedily combine non-overlapping hits and add resulting edge to graph
    """
    nolaps = None  # Coordinates for non-overlapping hits
    for aln_info in hit_stack:
        if not nolaps:
            nlp_bit = hit_stack[0][0]
            nlp_id_sum = hit_stack[0][1]
            nolaps = [hit_stack[0][2:]]
        else:            
            bit, id_sum, qbeg1, qend1, rbeg1, rend1 = aln_info  # Candidate hit
            for nlp_info in nolaps:
                qbeg2, qend2, rbeg2, rend2 = nlp_info  # Accepted hit
                if qbeg1 > qend2 or qend1 < qbeg2:  # No overlap in query
                    if rbeg1 > rend2 or rend1 < rbeg2:  # No overlap in ref
                        nlp_bit += bit
                        nlp_id_sum += id_sum

    # Make sure coverage (% identical bases) exceeds user-defined threshold
    nlp_cov = nlp_id_sum / min_len
    if nlp_cov >= min_cov:

        # G[A][B] != G[B][A] in nx Graphs, so add edges in node-sorted order
        id1, id2 = sorted([qry_id, ref_id])
        try:
            if nlp_bit > graph[id1][id2]['bit']:
                graph[id1][id2]['bit'] = nlp_bit
        except KeyError:
            graph.add_edge(id1, id2, bit=nlp_bit)


def add_top_hit_edges(graph, blast_handle, bscol=11, min_cov=0.5,
              allowed=False, disallowed=False):
    """Construct a graph from the top BLAST hits
    """
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

        if allowed:
            if not (graph.has_node(qry_id) and graph.has_node(ref_id)):
                continue

        elif disallowed:
            if qry_id in disallowed or ref_id in disallowed:
                continue

        bit = float(temp[bscol])
        id_frac = float(temp[2])/100
        aln_len = int(temp[3])
        id_sum = id_frac*aln_len
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
    """Compute average scores between and within each pair of organisms

    :param graph: Dictionary containing BLAST graph weighted with bit scores
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
            avgs.add_edge(orgA, orgB, cnt=1, sum=edata['bit'])

    for id1, id2, edata in avgs.edges(data=True):
        avgs[orgA][orgB]['avg'] = edata['sum']/edata['cnt']

    return avgs


def normalize_graph(graph, avgs, idchar):
    """Normalize graph using inter- & intra-organism averages

    This function inflates edge weights between less similar organisms, and
    deflates edge weights between more similar organism, converting all edge
    weights into dimensionless multiples of the overall average graph edge
    weight.

    :param graph: Dictionary containing BLAST graph weighted with bit scores
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
    """
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


if __name__ == "__main__":
    sys.exit(main())
