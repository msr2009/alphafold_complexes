"""
complexes_for_alphafold.py

a script that takes a single or multiple genes, isoforms, keywords, or fasta
files and outputs a alphafold-ready fasta file containing all combinations

Matt Rich, 12/2023
"""

import re
from itertools import groupby, combinations, product

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
		header.append(x.split()[0])
	header_str=">{}\n".format(":".join(header))
	seq_str = ":".join(zst[1])
	return header_str+seq_str

def main(protein_fasta, inputs, sizes):
	
	#first, read in fasta file containing all proteins
	proteins = {s[0]: s[1] for s in fasta_iter(protein_fasta)}
	
	#gather our sequences
	input_seqs = []
	
	#first if there are fasta:
	if inputs["FASTAS"] != None:
		for f in inputs["FASTAS"]:
			for s in fasta_iter(f):
				input_seqs.append([s])

	#we need to make regexs specifically for each input (except fasta)
	input_res = []
	#genes = "=[GENE] "
	if inputs["GENES"] != None:
		for x in inputs["GENES"]:
			input_res.append([re.compile("={} ".format(x)), None])
	#isoforms = "^ISOFORM "
	if inputs["ISOFORMS"] != None:
		for x in inputs["ISOFORMS"]:
			#check for ranges here
			if "[" in x:
				#split out range and put it the second element of these lists
				xr = x.split("[")[1][:-1].split("-")
				xres = x.split("[")[0]
				input_res.append([re.compile("^{} ".format(xres)), [int(xr[0])-1, int(xr[1])]])
			else:	
				input_res.append([re.compile("^{} ".format(x)), None])
	#keyword = "[KEYWORD]"
	if inputs["KEYWORDS"] != None:
		for x in inputs["KEYWORDS"]:
			input_res.append([re.compile("{}".format(x)), None])
	
#	#adjust sizes for max number
#	size_range = []
#	if sizes[0] == None:
#		size_range = sizes[2]
#	elif sizes[1] != None:
#		size_range = range(sizes[0], sizes[1]+1)
#	elif sizes[1] == None:
#		size_range = range(sizes[0], len(input_seqs))

	#search through proteins fasta to find matches
	
	for r in input_res:
		tmp_input = []
		for p in proteins:
			if r[0].search(p) != None:
				if r[1] == None:
					tmp_input.append((p, proteins[p]))
				else:
					tmp_input.append((p, proteins[p][r[1][0]:r[1][1]]))
		input_seqs.append(tmp_input)
	
	#now we have all our input sequences, we want to make 
	#all the relevant combinations of the proteins, choosing
	#up to one from each RE "bucket"
	#e.g., testing protein A + all isoforms of protein B individually
	

	#product(*x) gives us the combinations
	for c in product(*input_seqs):
		print(write_AF_fasta(c))



	
if __name__ == "__main__":
	
	from argparse import ArgumentParser

	parser = ArgumentParser()
	parser.add_argument('-p', '--proteins', action = 'store', type = str, dest = 'PROTEINS', 
		help = "fasta file containing all proteins to search")
	parser.add_argument('-g', '--gene', action = 'append', type = str, dest = 'GENES',
		help = 'search for gene name. returns all isoforms (e.g., unc-44)')
	parser.add_argument('-i', '--isoform', action = 'append', type = str, dest = 'ISOFORMS',
		help = 'search for specific isoform (e.g., B0350.2f)')
	parser.add_argument('-k', '--keyword', action = 'append', type = str, dest = 'KEYWORDS',
		help = 'search for keyword in description')
	parser.add_argument('-f', '--fasta', action = 'append', type = str, dest = 'FASTAS',
		help = 'input protein sequences as additional FASTA. No searching necessary!')
	parser.add_argument('--min', action = 'store', type = int, dest = "MIN_SIZE",
		help = 'minimum oligomeric size for output', default=1)
	parser.add_argument('--max', action = 'store', type = int, dest = "MAX_SIZE",
		help = 'maximum oligomeric size for output')
	parser.add_argument('-n', action = 'store', type = str, dest = "ARB_SIZE",
		help = 'comma-delimited list of arbitrary oligomeric numbers')

	args = parser.parse_args()
	
	combo_ns = [args.MIN_SIZE, None, None]
	if args.MAX_SIZE != None:
		combo_ns = [args.MIN_SIZE, args.MAX_SIZE, None]
	if args.ARB_SIZE != None:
		combo_ns = [None, None, [int(x) for x in args.ARB_SIZE.split(",")]]

	main(args.PROTEINS, {"GENES": args.GENES, "ISOFORMS": args.ISOFORMS, "KEYWORDS": args.KEYWORDS, "FASTAS": args.FASTAS}, combo_ns)
