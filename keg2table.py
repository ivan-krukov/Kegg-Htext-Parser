#Exports .keg BRITE hierarchy (htext) into a parsable table format

import re
import sqlite3
import tablib
from argparse import ArgumentParser

ec_pattern = re.compile(r"(\d+\.){3}\d")


"""KEGG htext stores level of hierarchy as a capital letter, starting at A
	This returns a numerical value corresponding to the level of hierarchy
"""
def hlevel(htext_line):
	return ord(htext_line[:1])-65

def hchar(level):
	return chr(level+65)

#It is expected that openning and closing sequences (a and b) are different
def between (s, a, b):
	return s[s.find(a)+len(a) : s.find(b)]

#It is expected that openning and closing sequences (a and b) are different
def rbetween (s, a, b):
	return s[s.rfind(a)+len(a) : s.rfind(b)]

#Strip the outermost pair of tags
def strip_tags(string):
	if ("<" in string and ">" in string):
		tag = between(string,"<",">")
		return between(string,"<{}>".format(tag),"</{}>".format(tag))
	else:
		return string

"""Line contains an ID, rest of line is a description
"""
def parse_id_line (content):
	tokens = content.split()
	return tokens[0]," ".join(tokens[1:])

"""Line contains an ID, an optional comma separated list ended with a semicolon and a description
"""
def parse_description_line (content):
	tokens = content.split()
	name = tokens[0]

	#try to find a comma separated list in the content
	optional_list = []
	list_end_index = 1
	for t in tokens[1:]:
		if t.endswith(","):
			optional_list.append(t[:-1])
		elif t.endswith(";"):
			optional_list.append(t[:-1])
			list_end_index += len(optional_list) 
			break
		else:
			break

	alternative_names = " ".join(optional_list)
	description = " ".join(tokens[list_end_index:])

	return name, alternative_names, description


"""Line contains an EC and a description
"""
def parse_ec_line (content):
	content = strip_tags(content)
	tokens = content.split()
	ec = tokens[0]
	description = " ".join(tokens[1:])
	return ec, description

"""Line contains an ID, description and a list with accessor in square brackets i.e. [ACC: 1 2 3]
"""
def parse_kegg_line (content):
	tokens = content.split()
	name = tokens[0]
	list_content = rbetween(content,"[","]")
	list_accessor = list_content[:list_content.find(":")]
	list_items = list_content[list_content.find(":")+1:]
	description = between(content,name,list_content)[1:-2]
	return name,description,list_accessor,list_items


class Stack:
	def __init__(self):
		self.stack = list()
	def push (self, item):
		self.stack.append(item)
	def pop (self):
		last = self.stack[-1]
		self.stack = self.stack[:-1]
		return last
	def top (self):
		return self.stack[-1]
	def content (self):
		return reversed(self.stack)
	def height (self):
		return len(self.stack)

class KeggParser:
	map_header = ["Gene name","Common names","Gene description","KO number","KO description","Map","KO mapping"]
	enzyme_header = ["4EC","4EC description","3EC","3EC description","2EC","2EC description","Enzyme type","Enzyme type description"]
	def __init__(self, encoding_line):
		self.kegg_map=False
		if encoding_line.startswith("+"):
			tokens = encoding_line.split()
			self.data_level = hlevel(tokens[0][1:])
			self.data_columns = tokens[1:]
			if tokens[1] == "Enzyme":
				self.enzyme_classification = True
			if len(self.data_columns) == 2:
				self.kegg_map = True

"""Some examples"""
a = "64144 Mllt1, AA407901, BAM11, ENL, LTG19; myeloid/lymphoid or mixed-lineage leukemia (trithorax homolog, Drosophila); translocated to, 1"
b = "17355 AF4/FMR2 family, member 1"
c = "214162 Mll1, 6430520K01, ALL-1, All1, Cxxc7, HRX, HTRX1, KIAA4050, KMT2A, Mll, mKIAA4050; myeloid/lymphoid or mixed-lineage leukemia 1 (EC:2.1.1.43)"
d = "K00567 methylated-DNA-[protein]-cysteine S-methyltransferase [EC:2.1.1.63]"

if __name__=="__main__":
	parser = ArgumentParser("Parses a '.keg' BRITE hierarchy file (htext)")
	parser.add_argument("keg_file",type=open,help="Input .keg file")
	parser.add_argument("output_file",type=str,help="Name of the output file or 'stdin'")
	parser.add_argument("-output_format",type=str,default="csv",
		choices=["xls","json","yaml","html","tsv","csv","sqlite"],help="Format of the output table (default .csv)")
	args = parser.parse_args()
	stack = Stack()
	header_line = list()

	for i, line in enumerate(args.keg_file):

		#first line is a data format header
		if (i == 0):
			parser = KeggParser(line)
			if parser.kegg_map:
				for column in parser.map_header:
					header_line.append(column)
				if parser.enzyme_classification:
					for column in parser.enzyme_header:
						header_line.append(column)
			else:
				for column in parser.data_columns:
					header_line.append(column)
				#python-fu: reverse traverse list without last element
				for category in range(parser.data_level-1,-1,-1):
					header_line.append(hchar(category))
			data = tablib.Dataset(headers=header_line)
			#print data.headers,len(data.headers)
			continue

		
		if (line.startswith("#") or line.startswith("%")):
			continue

		level = hlevel(line)

		#If we encounter a proper htext line:
		if level >= 0:
			content = line[1:].strip()
			#Make sure the line is not empty
			if content:
				table_line = list()
				#Close any previously opened categories as necessary
				while stack.height() and level <= hlevel(stack.top()):
					stack.pop()
				if level == parser.data_level:
					if (len(parser.data_columns) == 1):
						table_line.append(content)
					elif (parser.kegg_map):

						entity, classification = content.split("\t")
						for column in parse_description_line(entity):
							table_line.append(column) 
						for column in parse_kegg_line(classification):
							table_line.append(column)
					else:
						for field in content.split("\t"):
							table_line.append(field)
					for category in stack.content():
						if parser.enzyme_classification:
							for token in parse_ec_line(category[1:].strip()):
								table_line.append(token)
						else:
							table_line.append(category[1:].strip())
					#print table_line,len(table_line)
					data.append(table_line)
				else:
					stack.push(line)

	#TODO:handle STDIN as the output file
	open ("simple.xls","wb").write(data.xls)
	open(args.output_file,"wb").write(data.html)
			
