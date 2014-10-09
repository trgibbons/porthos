#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: tgibbons (tgibbons@umd.edu)
"""

# Standard Python libraries
import sys
import argparse


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

    graph = build_graph(blast_handle=args.blast, bscol=args.bscol-1)

    avgs = compute_organism_averages(graph=graph, idchar=args.idchar)

    if args.norm:
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

    parser.add_argument('blast',
                        nargs='?',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        help="Tab-delimited BLAST file (comment lines are " +
                             "okay) [def=stdin]")

    parser.add_argument('abc',
                        nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="'abc' output graph file [def=stdout]")

    args = parser.parse_args()

    return args


def build_graph(blast_handle, bscol=11):
    """Construct a graph from the top BLAST hits

    Bit scores from non-overlapping hits can be combined using simple addition,
    however the identification of the set of non-overlapping hits that provides
    the maximum score is an NP-complete problem. For this simple exemplary
    program in which I'm trying to avoid using non-standard Python modules, I
    decided not to tackle this aspect and am just using the top hits.

    In order to avoid using either a NetworkX graph or a pair of dictionaries,
    the query and subject identifiers are stored in lexicographical order.

    :param blast_handle: An open file handle containing non-self-alignments
        (can contain other alignments and/or comment lines beginning with a
        hash '#' character)
    :return graph: Dictionary containing BLAST graph weighted with bit scores
    """
    graph = dict()

    for line in blast_handle:
        temp = line.strip().split()
        if not temp:  # skip blank lines
            continue
        elif temp[0][0] == "#":  # skip comment lines
            continue
        elif temp[0] == temp[1]:  # skip self-hits
            continue

        id1, id2 = sorted([str(temp[0]), str(temp[1])])
        bit = float(temp[bscol])

        try:
            if bit > graph[id1][id2]:
                graph[id1][id2] = bit
        except KeyError:
            try:
                graph[id1][id2] = bit
            except KeyError:
                graph[id1] = dict(id2=bit)

    return graph


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

    for id1, nbrs in graph.iteritems():
        org1 = id1.split(idchar)[0]
        for id2, bit in nbrs.iteritems():
            org2 = id2.split(idchar)[0]
            orgA, orgB = sorted([org1, org2])

            try:
                avgs[orgA][orgB]['cnt'] += 1
                avgs[orgA][orgB]['sum'] += bit
            except KeyError:
                try:
                    avgs[orgA][orgB] = dict(cnt=1, sum=bit)
                except KeyError:
                    avgs[orgA] = dict(orgB=dict(cnt=1, sum=bit))

    for orgA, nbrs in avgs.iteritems():
        for orgB, stats in nbrs.iteritems():
            avgs[orgA][orgB]['avg'] = stats['sum']/stats['cnt']

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
    for id1, nbrs in graph.iteritems():
        org1 = id1.split(idchar)[0]
        for id2, bit in nbrs.iteritems():
            org2 = id2.split(idchar)[0]
            orgA, orgB = sorted([org1, org2])
            avg_bit = avgs[orgA][orgB]['avg']
            graph[id1][id2] /= avg_bit


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
    for id1, nbrs in graph.iteritems():
        for id2, bit in nbrs.iteritems():
            out_line = '\t'.join([id1, id2, str(bit)])+'\n'
            abc.write(out_line)


if __name__ == "__main__":
    sys.exit(main())
