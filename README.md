import Bio
import Bio.Seq

from Bio import SeqIO
from Bio.Seq import Seq

#For the path of the mutlifasta file

#filePATH = input("Please input a File Path, you may need to double the backslashes in order for program to work (ex: '\'--> '\\'): ")
WT_DNASeq_Input = input("Please provide an uppercase WT DNA sequence to use as a WT control: ")
WT_name = input("Please enter the name of your control protein: ")


#THE ASSUMPTIONS MADE FOR THIS SCRIPT TO WORK:
# - You receive great sequencing quality scores
# - Your sequence file is a Fasta or MultiFasta Format
# - Your sequence is less than ~1000 bps long
# - Your Identifiable region for all sequences is pETCON3/pET29 cutsight for NdeI or XhoI

# filePATH = ***YOUR FILE PATH HERE OR YOU CAN USE THE USER INPUT VARIABLE***

WT_DNASeq = WT_DNASeq_Input.upper()
bps = len(WT_DNASeq)

DNA_SEQ_OBJECT = Seq(WT_DNASeq)
WildType_AASeq = str(DNA_SEQ_OBJECT.translate())

# NdeI cutsight Identifier, Forward(Beginning of Sequence)
IdentifiersFW = "CATATG"
# XhoI cutsight Identifier, Reverse(End of Sequence)
IdentifiersRV = "CTCGAG"



# Should creat a list of tuples that contains the Seq ID and Seq Object(aka the sequence)
def List_of_Tupled_ID_Seq(filePATH):
    List_Of_Tuples = []
    MultiFastaFile = SeqIO.parse(filePATH, "fasta")
    for seq_record in MultiFastaFile:
        ActualID = seq_record.id
        ActualBPs = seq_record.seq
        Tupled_Seq_Info = ActualID, ActualBPs
        List_Of_Tuples.append(Tupled_Seq_Info)
    
    return List_Of_Tuples

# This should isolate the fragment of DNA starting at the Start Codon("ATG") after the Glycine linker and going until the length of the actual protein sequence

def BeginAlignment(IDs_Seqs):
    
    List_Of_FileIDs_Aligned_Sequences = []
    for seq in IDs_Seqs:
        if IdentifiersFW in seq[1]:
            a = seq[1].index(IdentifiersFW)
            List_Of_FileIDs_Aligned_Sequences.append((seq[0], str(seq[1][a+3:a+3+int(bps)].upper())))
        elif IdentifiersRV in seq[1]:
            a = seq[1].index(IdentifiersRV)
            List_Of_FileIDs_Aligned_Sequences.append((seq[0], str(seq[1][a-int(bps):a].upper())))  
        else:
            List_Of_FileIDs_Aligned_Sequences.append((seq[0],"TAATAATAA"))
            
    return List_Of_FileIDs_Aligned_Sequences

# Testing the output for the BeginAlignment, List_of_Tupled_ID_Seq functions that were created

Tupled_String_Info_Of_Sequences = BeginAlignment(List_of_Tupled_ID_Seq(filePATH))

# This function should output from the previous 2 functions aka take a the tupled info which contains start codon position 
# and then the basepairs that contain the protein sequence and translate them into Amino Acid Sequence with a plate ID
def Translate_String_Sequence(Tupled_String_SEQ):
    AminoAcid_Sequence_List_w_ID = []
    for item in Tupled_String_SEQ:
        FileID = item[0]
        SeqObject = Seq(item[1])
        AminoAcid_Sequence_List_w_ID.append((FileID, str(SeqObject.translate())))
    return AminoAcid_Sequence_List_w_ID


Variable_For_Translation_List = Translate_String_Sequence(Tupled_String_Info_Of_Sequences)
for item in Variable_For_Translation_List:
    print(len(item[1]))

def AA_Comparison_Operator(Amino_Acid_SEQ_w_ID):
# Part 1: Isolate variant Amino Acid sequences to later compare
    Only_Amino_Acid_Seqs = []
    for item in Amino_Acid_SEQ_w_ID:
        Only_Amino_Acid_Seqs.append((item[0], item[1]))

# Part 2: Record original Amino Acid, the Position and then Mutation        
    Record_Of_Mutations = []
    for AAseq in Only_Amino_Acid_Seqs:
        if AAseq[1] == WildType_AASeq:
            Record_Of_Mutations.append((AAseq[0],WT_name))
        elif AAseq[1] == "***":
            Record_Of_Mutations.append((AAseq[0], "Unable to Find Alignnment Identifier(s)"))
        else:
            try:
                for i in range(len(AAseq[1])):
                    if AAseq[1][i] != WildType_AASeq[i]:
                        Record_Of_Mutations.append((AAseq[0], WildType_AASeq[i], i+1 ,AAseq[1][i]))
            except Exception as e:
                Record_Of_Mutations.append((AAseq[0], e))

    return Record_Of_Mutations
            
    
AA_Comparison_Operator(Variable_For_Translation_List)

