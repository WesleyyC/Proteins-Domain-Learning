#!/usr/bin/python

import sys

# build the sequence dicn according to the BLOSUM matrix order
SEQ_2_NUM_DIC = {	'.': 0,
					'A': 1,
					'R': 2,
					'N': 3,
					'D': 4,
					'C': 5,
					'E': 6,
					'Q': 7,
					'G': 8,
					'H': 9,
					'O': 10,
					'I': 11,
					'L': 12,
					'K': 13,
					'M': 14,
					'F': 15,
					'P': 16,
					'U': 17,
					'S': 18,
					'T': 19,
					'W': 20,
					'Y': 21,
					'V': 22}

def get_New_File_Name():
	dash_Position = file_Name.find('_')
	colon_Position = file_Name.find(':')
	new_File_Name = file_Name[:dash_Position ]+"_toNum"+ file_Name[dash_Position:colon_Position+1]+ str(star_seq) + "-" + str(end_seq)+".csv"
	return new_File_Name

def get_Sequence():
	global seq_list, end_seq
	seq_list=[]

	f = open(file_Name,"r")

	start = False

	for line in f:
		for char in line:
			if not (char == ' ' or char == '\n'):
				seq_list.append(SEQ_2_NUM_DIC[char])

	if len(seq_list) < end_seq:
		end_seq = len(seq_list)

	f.close()

def write_New_Seq_File():
	new_File = open(get_New_File_Name(),"w")
	for num in seq_list[star_seq-1:end_seq]:
		new_File.write(str(num)+"\n")

	new_File.close()


if __name__ == "__main__":

	global file_Name, star_seq, end_seq
	file_Name = sys.argv[1]
	star_seq = int(sys.argv[2])
	end_seq = int(sys.argv[3])

	get_Sequence()

	write_New_Seq_File()
