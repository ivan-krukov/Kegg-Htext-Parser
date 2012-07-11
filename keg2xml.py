#Exports .keg BRITE hierarchy (htext) into a parsable xml format

from argparse import ArgumentParser
from cgi import escape

top = '<?xml version="1.0" encoding="UTF-8"?>\n<KEGG>'
bottom = "</KEGG>"

class Stack:
	def __init__(self):
		self.stack = []
	def push (self, item):
		self.stack.append(item)
	def pop (self):
		last = self.stack[-1]
		self.stack = self.stack[:-1]
		return last
	def top (self):
		return self.stack[-1]
	def height (self):
		return (len(self.stack))

"""KEGG htext stores level of hierarchy as a capital letter, starting at A
	This returns a numerical value corresponding to the level of hierarchy
"""
def hlevel(htext_line):
	return ord(htext_line[:1])-65

def hchar(level):
	return chr(level+65)

def print_with_header(string, header):
	print "<{}>".format(header)
	print string
	print "</{}>".format(header)

#It is expected that openning and closing sequences (a and b) are different
def between (s, a, b):
	return s[s.find(a)+len(a) : s.find(b)]

#It is expected that openning and closing sequences (a and b) are different
def rbetween (s, a, b):
	return s[s.rfind(a)+len(a) : s.rfind(b)]

#Strip the outermost pair of tags
def strip_tags(string):
	tag = between(string,"<",">")
	return between(string,"<{}>".format(tag),"</{}>".format(tag))

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



"""Some examples"""
a = "64144 Mllt1, AA407901, BAM11, ENL, LTG19; myeloid/lymphoid or mixed-lineage leukemia (trithorax homolog, Drosophila); translocated to, 1"
b = "17355 AF4/FMR2 family, member 1"
c = "214162 Mll1, 6430520K01, ALL-1, All1, Cxxc7, HRX, HTRX1, KIAA4050, KMT2A, Mll, mKIAA4050; myeloid/lymphoid or mixed-lineage leukemia 1 (EC:2.1.1.43)"
d = "K00567 methylated-DNA-[protein]-cysteine S-methyltransferase [EC:2.1.1.63]"

if __name__=="__main__":
	parser = ArgumentParser("Parses a .keg BRITE hierarchy file (htext)")
	parser.add_argument("keg_file",type=str,help="Input .keg file")
	args = parser.parse_args()

	stack = Stack()
	
	print top
	for line in open(args.keg_file):
		level = hlevel(line)

		#TODO: resolve comments and other things - need to look up more htext files
		if level >= 0:
			content = escape(line[1:].strip())
			if content:
				while stack.height() and stack.top() >= level:
					print "</"+hchar(stack.pop())+">"	
				stack.push(level)
				print "<"+hchar(level)+">"
				print content

	#unwind the stack
	while stack.height() > 0:
		print "</"+hchar(stack.pop())+">"
	print bottom
					
