#!/usr/bin/python

import re
import sys
import os
import getopt


def main():
	params = parseArgs()
	
	if params.fasta:
		header = False
		headerline=""
		num=1
		out=(os.path.splitext(params.fasta)[0])+".tsv"
		with open(out, 'w') as OFH:
			for samp in read_fasta(params.fasta):
				info = samp[0].split("_")
				seq = samp[1]
				
				name = "".join(info[:-1])
				if params.append:
					name = name + "_" + str(num)
					num+=1
				split_dat = info[-1].replace("["," ").replace("]"," ").split(" ")
				
				oline = name + "\t" + split_dat[0] + "\t"
				if len(split_dat) > 1:
					oline = oline + split_dat[1]+ "\t"
					oline = oline + split_dat[3]+"\t"
				oline = oline + seq + "\n"
						
				
				if header == False:
					headerline = 'name\taccession\t'
					if len(split_dat) > 1:
						headerline = headerline + "lat\tlong\t"
					headerline = headerline+"seq\n"
					header=True
					OFH.write(str(headerline))
				
				OFH.write(str(oline))
		OFH.close()
	elif params.table:
		header = False
		headerline=""
		coords=False
		out=(os.path.splitext(params.table)[0])+".fasta"
		with open(out, 'w') as OFH:
			for line in read_table(params.table):
				if not header:
					header=True
					if len(line) > 3:
						coords=True
					continue
				fas_header = ">" + line[0] + "_" + line[1]
				if coords:
					fas_header = fas_header + "[" + str(line[2]) + " N "
					fas_header = fas_header + str(line[3]) + " W]"
				fas_header = fas_header + "\n"
				OFH.write(fas_header)
				
				seq = line[-1] + "\n"
				OFH.write(seq)
		OFH.close()
				
				
def read_table(fas):
	if os.path.exists(fas):
		with open(fas, 'r') as fh:
			try:
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					yield(line.split("\t"))
			except IOError:
				print("Could not read file ",fas)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%fas)
		
#Read samples as FASTA. Generator function
def read_fasta(fas):

	if os.path.exists(fas):
		with open(fas, 'r') as fh:
			try:
				contig = ""
				seq = ""
				for line in fh:
					line = line.strip()
					if not line:
						continue
					#print(line)
					if line[0] == ">": #Found a header line
						#If we already loaded a contig, yield that contig and
						#start loading a new one
						if contig:
							yield([contig,seq]) #yield
							contig = "" #reset contig and seq
							seq = ""
						contig = (line.replace(">",""))
					else:
						seq += line
				#Iyield last sequence, if it has both a header and sequence
				if contig and seq:
					yield([contig,seq])
			except IOError:
				print("Could not read file ",fas)
				sys.exit(1)
			finally:
				fh.close()
	else:
		raise FileNotFoundError("File %s not found!"%fas)

#Object to parse command-line arguments
class parseArgs():
	def __init__(self):
		#Define options
		try:
			options, remainder = getopt.getopt(sys.argv[1:], 'nf:t:h', \
			["help"])
		except getopt.GetoptError as err:
			print(err)
			self.display_help("\nExiting because getopt returned non-zero exit status.")
		#Default values for params
		#Input params
		self.fasta=None
		self.table=None
		self.append=False


		#First pass to see if help menu was called
		for o, a in options:
			if o in ("-h", "-help", "--help"):
				self.display_help("Exiting because help menu was called.")

		#Second pass to set all args.
		for opt, arg in options:
			arg = arg.strip()
			opt = opt.replace("-","")
			#print(opt,arg)
			if opt == 'f':
				self.fasta = arg
			elif opt == "t":
				self.table=arg
			elif opt=="n":
				self.append=True
			elif opt == 'h' or opt == 'help':
				pass
			else:
				assert False, "Unhandled option %r"%opt

		#Any changes necessary
		if not self.fasta and not self.table:
			self.display_help("must input one of <-t> or <-f>")


	def display_help(self, message=None):
		if message is not None:
			print()
			print (message)
		print ("\nfetcher.py\n")
		print("Author: Tyler K Chafin, University of Arkansas")
		print ("Contact: tkchafin@uark.edu")
		print ("Description: Convert autoStreamTree fasta to a table and back")
		print("""
	Options:
		-f	: Input fasta to convert to table
		-t	: Input table to convert to fasta
		-n	: Append unique numbering to each sample ID
		-h,--help	: Displays help menu

""")
		print()
		sys.exit()

#Call main function
if __name__ == '__main__':
    main()
