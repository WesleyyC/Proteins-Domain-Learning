#!/usr/bin/python

import sys
from random import randint

def write_Random_Seq_File(seq_length):
	new = open("rand_seq.csv",'w')
	id = 0
	seq = 0
	count = 0
	x = 0
	y = 0
	z = 0
	fill_id=1;

	for id in range(seq_length)
		new.write(str(id)+","+str(randint(1,20)))
	
	new.close()


if __name__ == "__main__":

	seq_length= sys.argv[1]

	write_New_Seq_File(seq_length)




