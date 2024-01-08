#PACKAGES USED
import Bio
import Bio.Seq

from Bio import SeqIO
from Bio.Seq import Seq

#THE ASSUMPTIONS MADE FOR THIS SCRIPT TO WORK:
# - You receive great sequencing quality scores
# - Your sequence file is a Fasta or MultiFasta Format
# - Your sequence is less than ~1000 bps long
# - Your Identifiable region for all sequences is pETCON3/pET29 cutsight for NdeI or XhoI

# Part 1: Read multiFASTA File
# Part 2: Isolate variant Amino Acid sequences to later compare
# Part 3: Record original Amino Acid, the Position and then Mutation        


