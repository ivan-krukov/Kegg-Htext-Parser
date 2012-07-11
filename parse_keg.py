#Parses the KEG htext file format

from argparse import ArgumentParser
from collections import defaultdict
from collections import namedtuple
import re

"""Hierarchy structure:
	A single Entry corresponds to an EC number. 
		EC number has a description and a list of genes.
	Each gene has a unique sequence name, optional common name and a description.
		To each protein, a K number (KEGG ID) is mapped
	A KEGG ID consists of a name and a list of EC numbers that reference that protein
"""

Entry = namedtuple("Entry", ["ec", "description", "genes"])
Gene = namedtuple("Gene", ["sequence_name", "common_name", "description", "kegg_id"])
Kegg_ID = namedtuple("Kegg_ID", ["name", "description", "ec_numbers"])

def between(string, char_a, char_b):
	return string [string.find(char_a)+1 : string.find(char_b)]

if __name__=="__main__":

	parser = ArgumentParser("Parses a .keg BRITE hierarchy file (htext)")
	parser.add_argument("keg_file",type=str,help="Input .keg file")
	args = parser.parse_args()

	ec_pattern = re.compile(r"([-0-9]+\.){3}[-0-9]")

	hierarchy = defaultdict(Entry)
	line_buffer = []
	ec_line = ""

	for line in open(args.keg_file):
		
		if line.startswith("D") and not line_buffer:	
			ec_line = line

		elif line.startswith("E") and ec_line:
			line_buffer.append(line)

		elif line.startswith("D") and line_buffer:
			ec_line_tokens = ec_line.split()
			ec_number = ec_line_tokens[1]
			ec_description = " ".join(ec_line_tokens[2:])
				
			print ec_number,ec_description
			for description_line in line_buffer:
				protein_entry,kegg_entry = description_line.split("\t")

				protein_entry_tokens = protein_entry.split()
				protein_sequence_name = protein_entry_tokens[1]
				if (";" in protein_entry_tokens[2]):
					protein_common_name = protein_entry_tokens[2][:-1]
					protein_description = " ".join(protein_entry_tokens[3:])
				else:
					protein_common_name = ""
					protein_description = " ".join(protein_entry_tokens[2:])
				
				kegg_entry,kegg_mapped_ecs = kegg_entry.split("[")
				kegg_entry_tokens = kegg_entry.split()
				kegg_id = kegg_entry_tokens[0]
				kegg_description = kegg_entry_tokens[1:]

				kegg_ecs = ec_pattern.findall(kegg_mapped_ecs)
				print """\tID: {seq_name}; Common: {name}\n\t{description}\n\tKEGG: {K_id}""".format(
					seq_name = protein_sequence_name, 
					name = protein_common_name, 
					description = protein_description, 
					K_id = kegg_id)

			ec_line = line
			line_buffer = []
			
