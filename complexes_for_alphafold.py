"""
complexes_for_alphafold.py

a script that takes a single or multiple genes, isoforms, keywords, or fasta
files and outputs a alphafold-ready fasta file containing all combinations

Matt Rich, 12/2023
"""

import re
from itertools import groupby, combinations, product

mandatory_string = "%"

def fasta_iter(fasta_name):
	"""
	modified from Brent Pedersen
	https://www.biostars.org/p/710/
	given a fasta file. yield tuples of header, sequence
	"""
	# first open the file outside 
	fh = open(fasta_name)

	# ditch the boolean (x[0]) and just keep the header or sequence since	
	# we know they alternate.
	faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))

	for header in faiter:
		# drop the ">"
		headerStr = header.__next__()[1:].strip()
	
		# join all sequence lines to one.
		seq = "".join(s.strip() for s in faiter.__next__())

		yield (headerStr, seq)

def write_AF_fasta(seq_tuples):
	header = []
	seq = []
	#zip the seq tuples into header and seq tuples
	zst = list(zip(*seq_tuples))
	for x in zst[0]:
		header.append(x.split()[0].strip(","))
	header_str=">{}\n".format(":".join(header))
	seq_str = ":".join(zst[1])
	return header_str+seq_str

def parse_complex_sizes(s):
	#s is comma-delimited, either an int or a range
	items = s.split(",")
	ranges = []
	for item in items:
		if "-" in item:
			start, end = map(int, item.split("-"))
			ranges.extend(range(start, end+1))
		else:
			ranges.append(int(item))
	return iter([x-1 for x in ranges])

def main(protein_fasta, inputs, sizes):
	
	#first, read in fasta file containing all proteins
	proteins = {s[0]: s[1] for s in fasta_iter(protein_fasta)}
	
	#gather our sequences
	input_seqs = []

	#keep track of which are mandatory in the output
	mandatories = []
	
	#first if there are fasta:
	if inputs["FASTAS"] != None:
		for f in inputs["FASTAS"]:
			for s in fasta_iter(f.rstrip(mandatory_string)):
				if f.endswith(mandatory_string):
					mandatories.append(s)
				input_seqs.append([s])

	#we need to make regexs specifically for each input (except fasta)
	input_res = []
	#genes = "=[GENE] "
	if inputs["GENES"] != None:
		for x in inputs["GENES"]:
			input_res.append([re.compile("={} ".format(x)), None, False])
	#isoforms = "^ISOFORM "
	if inputs["ISOFORMS"] != None:
		for x in inputs["ISOFORMS"]:
			mand = False
			if x.endswith(mandatory_string):
				mand = True
				x = x.rstrip(mandatory_string)
			#check for ranges here
			if "[" in x:
				#split out range and put it the second element of these lists
				xr = x.split("[")[1].split("]")[0].split("-")
				xres = x.split("[")[0]
				input_res.append([re.compile("^{} ".format(xres)), [int(xr[0])-1, int(xr[1])], mand])
			else:	
				input_res.append([re.compile("^{} ".format(x)), None, mand])
	#keyword = "[KEYWORD]"
	if inputs["KEYWORDS"] != None:
		for x in inputs["KEYWORDS"]:
			input_res.append([re.compile("{}".format(x)), None, False])
	
	#search through proteins fasta to find matches
	#for inputs 

	for r in input_res:
		tmp_input = []
		for p in proteins:
			if r[0].search(p) != None:
				if r[1] == None:
					tmp_input.append((p, proteins[p]))
				else:
					tmp_input.append((p, proteins[p][r[1][0]:r[1][1]]))
				#if this is a mandatory output
				if r[2]:
					mandatories.extend(tmp_input)
		input_seqs.append(tmp_input)
	
#	print("INPUTSEQS: {}".format(input_seqs))
#	print("MANDATORY: {}".format(mandatories))
	#now we have all our input sequences, we want to make 
	#all the relevant combinations of the proteins, choosing
	#up to one from each RE "bucket"
	#e.g., testing protein A + all isoforms of protein B individually
	#I want to enumerate all possible combinations from each bucket
	#then filter for what we really want.
	
	#if the user defines what complex sizes they want
	#then make that iterator
	complex_sizes = range(len(input_seqs))
	if sizes != None:
		complex_sizes = parse_complex_sizes(sizes)
	
	#then enumerate all possible with those sizes
	output_seqs = []
	for i in complex_sizes:
		for x in combinations(range(len(input_seqs)), i+1):
			for combo in product(*[ input_seqs[y] for y in x ]):
				if sum([ c in mandatories for c in combo ]) == len(mandatories):
					print(write_AF_fasta(combo))
	
if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-p', '--proteins', action = 'store', type = str, dest = 'PROTEINS', 
		help = "fasta file containing all proteins to search", required=True)
	parser.add_argument('-g', '--gene', action = 'append', type = str, dest = 'GENES',
		help = 'search for gene name. returns all isoforms (e.g., unc-44)')
	parser.add_argument('-i', '--isoform', action = 'append', type = str, dest = 'ISOFORMS',
		help = 'search for specific isoform (e.g., B0350.2f). Can specify \
				residues using [start-stop]. Append "{}" to make mandatory in complexes.'.format(mandatory_string))
	parser.add_argument('-k', '--keyword', action = 'append', type = str, dest = 'KEYWORDS',
		help = 'search for keyword in description')
	parser.add_argument('-f', '--fasta', action = 'append', type = str, dest = 'FASTAS',
		help = 'input protein sequences as FASTA.')
	parser.add_argument('-o', '--oligomer', action = 'store', type = str, dest = 'OLIGO',
		help = 'size of complexes to output. Can be range (1-4) or comma-delimited (1,3,5). Default contains one of each input.',
		default=None)

	args = parser.parse_args()
	
	main(args.PROTEINS, {"GENES": args.GENES, "ISOFORMS": args.ISOFORMS,
			"KEYWORDS": args.KEYWORDS, "FASTAS": args.FASTAS}, args.OLIGO)
