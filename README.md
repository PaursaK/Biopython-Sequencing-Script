import Bio
import Bio.Seq

from Bio import SeqIO
from Bio.Seq import Seq

#For the path of the mutlifasta file

filePATH = input("Please input a File Path, you need to double the backslashes in order for program to work (ex: '\\'--> '\\\'): ")
bps = input("Please provide the raw basepairs of the WT protein sequence you are investigating. This exclude the stop codon and yeast display components: ")
WildType_Seq = input("Please provide an uppercase amino acid sequence to use as a WT comparison: ")
Identifiers = input("Please provide uppercase basepair sequence of that should remains unchanged in all constructs looking to be compared: ")
WT_name = input("Please enter the name of your control protein: ")

#filePATH = "C:\\Users\\Paursa Kamalian\\Downloads\\DAVIDSHOULTZ10-17-2022_015359_651\\30-774315135.fasta"
#bps = 354
#WildType_Seq = "MSEEQIRQFLRRFYEALDSGDADTAASLFHPGVTIHLWDGVTFTSREEFREWFERLFSTSKDAQREIKSLEVRGDTVEVHVQLHATHNGQKHTVDLTHHWHFRGNRVTEVRVHINPTG"
#Identifiers = "ATGAGC"
#WT_name = "WT LUXSIT"




#Should creat a list of tuples that contains the Seq ID and Seq Object(aka the sequence)
def List_of_Tupled_ID_Seq(filePATH):
    List_Of_Tuples = []
    MultiFastaFile = SeqIO.parse(filePATH, "fasta")
    for seq_record in MultiFastaFile:
        ActualID = seq_record.id
        ActualBPs = seq_record.seq
        Tupled_Seq_Info = ActualID, ActualBPs
        List_Of_Tuples.append(Tupled_Seq_Info)
    
    return List_Of_Tuples



#Should isolate list of sequences from the tuple of ID and convert the Seq Object into a string
def ExtractDNA(Tupled_SEQList_Info):
    list_of_seqs = []
    for item in Tupled_SEQList_Info:
        list_of_seqs.append(str(item[1]).upper())

    return list_of_seqs

#This should isolate the fragment of DNA starting at the Start Codon("ATG") after the Glycine linker and going until the length of the actual protein sequence
#Need to make this more robust in order to be utilized for any type of protein we are trying to sequence that is less than 1000bps and we get sequencing results
#Right now this is tailored towards luxsit

def BeginAlignment(list_of_sequences):
    
    list_of_indexes = []
    for seq in list_of_sequences:
        if Identifiers in seq:
            a = seq.index(Identifiers)
            list_of_indexes.append((a, seq[a:a+int(bps)],len(seq[a:a+int(bps)])))
        elif "N" or "-" in seq:
            list_of_indexes.append("TAATAATAA")
            
    return list_of_indexes

#Testing the output for the BeginAlignment, ExtractDNA and List_of_Tupled_ID_Seq functions that were created

Tupled_String_Info_Of_Sequences = BeginAlignment(ExtractDNA(List_of_Tupled_ID_Seq(filePATH)))

#This function should output from the previous 3 functions aka take a the tupled info which contains start codon position 
#and then the basepairs that contain the protein sequence and translate them into Amino Acid Sequence with a plate ID
def Translate_String_Sequence(Tupled_String_SEQ):
    AminoAcid_Sequence_List_w_ID = []
    for item in Tupled_String_SEQ:
        SeqObject = Seq(item[1])
        AminoAcid_Sequence_List_w_ID.append((Tupled_String_SEQ.index(item)+1,str(SeqObject.translate())))
    return AminoAcid_Sequence_List_w_ID


Variable_For_Translation_List = Translate_String_Sequence(Tupled_String_Info_Of_Sequences)

#for L in Variable_For_Translation_List:
#    print(L)
    

def AA_Comparison_Operator(Amino_Acid_SEQ_w_ID):
#Part 1: Isolate variant Amino Acid sequences to later compare
    Only_Amino_Acid_Seqs = []
    for item in Amino_Acid_SEQ_w_ID:
        Only_Amino_Acid_Seqs.append(item[1])

#Part 2: Record original Amino Acid, the Position and then Mutation        
    Record_Of_Mutations = []
    idx = 0
    for AAseq in Only_Amino_Acid_Seqs:
        idx += 1
        if AAseq == WildType_Seq:
            Record_Of_Mutations.append((idx,WT_name))
        elif AAseq == "TAATAATAA":
            Record_Of_Mutations.append((idx, "Needs manual review or sequencing needs to be repeated"))
        else:
            for i in range(len(AAseq)):
                if AAseq[i] != WildType_Seq[i]:
                    Record_Of_Mutations.append((idx, WildType_Seq[i], i+1 ,AAseq[i]))

    return Record_Of_Mutations
            
    
AA_Comparison_Operator(Variable_For_Translation_List)

