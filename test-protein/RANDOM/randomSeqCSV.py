#!/usr/bin/python

import sys
from random import randint


def write_Random_Seq_File(seq_length):
    new = open("rand_seq_" + str(seq_length) + ".csv", 'w')

    for id in range(seq_length):
        new.write(str(id) + "," + str(randint(1, 20)) + "\n")

    new.close()


if __name__ == "__main__":

    seq_length = int(sys.argv[1])

    write_Random_Seq_File(seq_length)
