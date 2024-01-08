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



# Should creat a list of tuples that contains the Seq ID and Seq Object(aka the sequence)


# This should isolate the fragment of DNA starting at the Start Codon("ATG") after the Glycine linker and going until the length of the actual protein sequence


# Testing the output for the BeginAlignment, List_of_Tupled_ID_Seq functions that were created


# This function should output from the previous 2 functions aka take a the tupled info which contains start codon position 
# and then the basepairs that contain the protein sequence and translate them into Amino Acid Sequence with a plate ID

# Part 1: Isolate variant Amino Acid sequences to later compare


# Part 2: Record original Amino Acid, the Position and then Mutation        


