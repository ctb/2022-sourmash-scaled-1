#! /usr/bin/env python
"""Calculate Jaccard similarities between genomes at a scaled=1.

This code outputs a numpy matrix that is compatible with 'sourmash plot'.

CTB note to future self: we are using SourmashSignature rather than
MinHash because it's much, much faser to run
SourmashSignature.add_sequence! The underlying reason is that
signatures constructed this way use the BTree implementation of MinHash.
"""
import sourmash
import sys
import argparse
import screed
import numpy as np

from sourmash.command_sketch import _signatures_for_sketch_factory


def consume(sig, filename):
    print(f"consuming '{filename}'")
    for record in screed.open(filename):
        print('record:', record.name)
        sig.add_sequence(record.sequence, force=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('genomes', nargs='+')
    p.add_argument('-k', '--ksize', default=31, type=int)
    p.add_argument('-o', '--output-prefix',
                   help="output numpy matrix/labels for use with sourmash plot")
    args = p.parse_args()

    N = len(args.genomes)
    print(f"Got {len(args.genomes)} genome files.")

    param_str = f'k={args.ksize},scaled=1'
    print("using:", param_str)

    factory = _signatures_for_sketch_factory([param_str], 'dna')

    siglist = []
    for filename in args.genomes:
        # make new SourmashSignature
        ss = factory()[0]

        # consume!
        consume(ss, filename)

        # store
        siglist.append(ss)

    # build a numpy matrix
    mat = np.ones((N, N))

    # fill it in!
    for i in range(N):
        ss_i = siglist[i]
        for j in range(i):
            ss_j = siglist[j]
            mat[i, j] = ss_i.jaccard(ss_j)
            mat[j, i] = mat[i, j] # Jaccard is symmetric

        mat[i, i] = 1           # jaccard of self to self is always 1
        
    with open(args.output_prefix, "wb") as fp:
        np.save(fp, mat)

    with open(f"{args.output_prefix}.labels.txt", "wt") as fp:
        for g in args.genomes:
            print(g, file=fp)

    np.set_printoptions(precision=3, suppress=True)
    print(mat)


if __name__ == '__main__':
    sys.exit(main())
