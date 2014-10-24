#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 22:01:00 2014

@author: tgibbons
"""

# Standard Python libraries
from sys import stderr
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

    if args.allowed:
        allowed = set(line.strip() for line in open(args.allowed))
    else:
        allowed = False

    if args.disallowed:
        disallowed = set(line.strip() for line in open(args.disallowed))
    else:
        disallowed = False

    for edge in args.abc_in:
        temp = edge.strip().split()
        if not temp:  # reprint blank lines
            continue
        elif temp[0][0] == "#":  # reprint comment lines
            args.abc_out.write(edge)
            continue

        id1 = str(temp[0])
        id2 = str(temp[1])

        if allowed:
            if not (id1 in allowed and id2 in allowed):
                continue

        if disallowed:
            if id1 in disallowed or id2 in disallowed:
                continue

        args.abc_out.write(edge)


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
        description="Filter edges in an ABC-format graph file")

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

    parser.add_argument('abc_in',
                        nargs='?',
                        type=argparse.FileType('r'),
                        default=sys.stdin,
                        help="ABC-formatted input graph file [def=stdin]")

    parser.add_argument('abc_out',
                        nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help="ABC-formatted output graph file [def=stdout]")

    args = parser.parse_args()

    return args


if __name__ == "__main__":
    sys.exit(main())
