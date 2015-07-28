#!/usr/bin/python

import sys

AMINO_IDX = 3;
R_ID_IDX = 5

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

TRIP_2_SING_DIC = {	'ALA': 'A',
					'ARG': 'R',
					'ASN': 'N',
					'ASP': 'D',
					'CYS': 'C',
					'GLU': 'E',
					'GLN': 'Q',
					'GLY': 'G',
					'HIS': 'H',
					'HYP': 'O',
					'ILE': 'I',
					'LEU': 'L',
					'LYS': 'K',
					'MET': 'M',
					'PHE': 'F',
					'PRO': 'P',
					'GLP': 'U',
					'SER': 'S',
					'THR': 'T',
					'TRP': 'W',
					'TYR': 'Y',
					'VAL': 'V'}


def write_New_Seq_File():
	org = open(file_Name,"r")
	new = open("new_"+ file_Name,'w')
	id = 0
	seq = 0
	count = 0
	x = 0
	y = 0
	z = 0

	for line in org:
		if "ATOM" in line:
			new_id=int(line[22:26].strip())
			if (new_id == id):
				count = count + 1
				x = x + float(line[30:38].strip())
				y = y + float(line[38:46].strip())
				z = z + float(line[46:53].strip())
			else:
				if(count!=0):
					new.write(str(id)+","+str(seq)+","+str(x/count)+","+str(y/count)+","+str(z/count)+"\n")
				id = new_id
				seq = SEQ_2_NUM_DIC[TRIP_2_SING_DIC[line[17:20].strip()]]
				count = 1
				x = float(line[30:38].strip())
				y = float(line[38:46].strip())
				z = float(line[46:53].strip())

	new.write(str(id)+","+str(seq)+","+str(x/count)+","+str(y/count)+","+str(z/count)+"\n")
	org.close()
	new.close()


if __name__ == "__main__":

	global file_Name
	file_Name = sys.argv[1]

	write_New_Seq_File()




