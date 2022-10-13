# Biopython-Sequencing-Script

import Bio
import Bio.Seq

from Bio import SeqIO
from Bio.Seq import Seq



filePATHtest =  "C:\\Users\\Paursa Kamalian\\Downloads\\DAVIDSHOULTZ10-13-2022_054049_084\\30-774547905.fasta"

#Should creat a list of tuples that contains the Seq ID and Seq Object(aka the sequence)
def List_of_Tupled_ID_Seq(FilePATH):
    List_Of_Tuples = []

    for seq_record in SeqIO.parse(filePATH, "fasta"):
        ActualID = seq_record.id
        ActualBPs = seq_record.seq
        Tupled_Seq_Info = ActualID, ActualBPs
        List_Of_Tuples.append(Tupled_Seq_Info)
    
    return List_Of_Tuples     

#Should isolate list of sequences from the tuple of ID and Seq Objects
def ExtractDNA(Tupled_SEQList_Info):
  list_of_seqs = []
  for item in Tupled_Info:
    list_of_seqs.append(item[1])

  return list_of_seqs

#Creates a string of the sequence object from the fasta files
string_seq_list = []
for sequence in list_of_seqs:
    newseq = str(sequence)
    string_seq_list.append(newseq)
    
