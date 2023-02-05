# import required libraries
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import re
import sys

# Parse command line arguments
parser = ArgumentParser(formatter_class=ArgumentDefaultsHelpFormatter)
parser.add_argument("nargs", nargs="+", help="fastq files")
parser.add_argument("--quality_filter", default=45, type=int, help="filter quality score")
args = vars(parser.parse_args())

class FastqParser:
    '''This is a class defined to read a fastq file, 
    parse it and generate a faster file with the sequence 
    ID and the raw sequence'''
    
    # Set up parameters
    fastq_files = args["nargs"]
    quality_filter = args["quality_filter"]
    
    # create a dictionary for ascii encoding characters representing quality score for each base
    qscore_dict = {"!": 33, "\"": 34, "#": 35, "$": 36, "%": 37, "&": 38,
     "'": 39, "(": 40, ")": 41, "*": 42, "+": 43, ",": 44, "-": 45, 
     ".": 46, "/": 47, "0": 48, "1": 49, "2": 50, "3": 51, "4": 52, "5": 53, 
     "6": 54, "7": 55, "8": 56, "9": 57, ":": 58, ";": 59, "<": 60, "=": 61,
     ">": 62, "?": 63, "@": 64, "A": 65, "B": 66, "C": 67, "D": 68, "E": 69,
     "F": 70, "G": 71, "H": 72, "I": 73, "J": 74, "K": 75, "L": 76, "M": 77,
     "N": 78, "O": 79, "P": 80, "Q": 81, "R": 82, "S": 83, "T": 84, "U": 85,
     "V": 86, "W": 87, "X": 88, "Y": 89, "Z": 90, "[": 91, "\\": 92, "]": 93, 
     "^": 94, "_": 95, "`": 96, "a": 97, "b": 98, "c": 99, "d": 100, 
     "e": 101, "f": 102, "g": 103, "h": 104, "i": 105, "j": 106, "k": 107, 
     "l": 108, "m": 109, "n": 110, "o": 111, "p": 112, "q": 113, "r": 114,
     "s": 115, "t": 116, "u": 117, "v": 118, "w": 119,
     "x": 120, "y": 121, "z": 122, "{": 123, "|": 124, "}": 125, "~": 126}

    # create a constructor for the class FastqParser
    def __init__(self, fastq_file):
        self.fastq_file = fastq_file

    def check_files(self, fastq_file):
        '''this method is used to check the fastq files for bad entry'''
        self.fastq_file = fastq_file
        file = open(self.fastq_file) # open the fastq file and read line by line
        fastq_lines = file.readlines()
        file.close()
        # group the fastq file into a block of four lines for each sequence
        seqs = [fastq_lines[i:i + 4] for i in range(0, len(fastq_lines), 4)]
        self.fastq_file = seqs
        bases = ['A', 'C', 'T', 'G', 'N'] # a list of the possible nucleotides
        for seq in self.fastq_file:
            ID =re.findall(r"(?<=@)[A-Za-z0-9-\.]+", seq[0])[0] # regex pattern to match sequence ID
            # check that each fastq sequence exist in a group of four lines
            if len(seq) < 4:
                print(f'{fastq_file} {ID} Error: Fastq sequence does exist in a block of 4 lines')
            # check that each fastq sequence contains normal bases i.e. A,C,G,T and N for missing call
            error_bases = []
            for char in seq[1].strip():
                if char not in bases:
                    error_bases.append(char)
            if error_bases:
                print(f'{fastq_file} {ID} Error: Invalid bases in the sequence - {error_bases}')
        return self.fastq_file
    
    def get_identifier(self, seq):
        '''get the sequence identifier for each 
        fastq sequence in the source fastq file'''
        self.seq = seq
        ID =re.findall(r"(?<=@)[A-Za-z0-9-\.]+", self.seq[0])[0] # regex pattern to match sequence ID
        return ID

    def missing_base(self, seq):
        '''determine the count of missing bases
        in each fastq sequence'''
        self.seq = seq
        count = 0 # initialize count of missing base to 0
        missing_base = 'N'
        for char in self.seq[1].strip():
            if char == missing_base: # if a nucleotide base is N, then count is incremented by 1
                count += 1
        return count
    
    def percent_GC(self, seq):
        '''determine the percentage GC of each sequence'''
        self.seq = seq
        C, G, A, T = 0, 0, 0, 0 # initialize the count of C,G,A,T to 0
        for nuc in self.seq[1]: # for each nucleotide in the sequence string 
            if nuc == 'C':
                C += 1
            elif nuc == 'G':
                G += 1
            elif nuc == 'T':
                T += 1  
            elif nuc == 'A':
                A += 1
        sum_of_nuc = C + G + T + A # add the count of each bases in the sequence together
        GC = G + C # add count of G and C
        per_GC = (GC / sum_of_nuc) * 100 # obtain percentage GC
        return per_GC
    
    def q_score(self, seq):
        '''determine the average quality score of each fastq sequence. 
        This is necessary for quality filtering'''
        self.seq = seq
        scores = [] # create an empty list to store the score of each nucleotide base of the sequence
        for char in self.seq[3].strip(): # loop through the quality symbol encoding quality score
            score = FastqParser.qscore_dict[char] # get the ascii score of the quality symbol
            scores.append(score) # append the score to the list of scores for the sequence
        qscore = sum(scores)/len(scores) # obtain the average quality score
        return qscore

# the code below will parse input fastq files based on the class module
for file in FastqParser.fastq_files: # loop through the fastq files provided in the command line
    filepath = file[:-6] + '.fasta' # generate the same name for the output fasta file
    fastq_parser = FastqParser(file) # create an object of the class "FastqParser"
    seqs = fastq_parser.check_files(file) # group the fatsq file into block of four lines (reads)
    
    for seq in seqs: # loop through each sequence (each block of four lines containing each fastq read)
        ID = fastq_parser.get_identifier(seq) # generate the sequence ID for the read
        missing_count = fastq_parser.missing_base(seq) # determine the number of missing call for the read
        per_GC = fastq_parser.percent_GC(seq) # obtain the percentage GC
        qscore = fastq_parser.q_score(seq) # obtain the average quality score of the read
        
        # print the filename, sequence ID, missing base count and %GC on the command line
        print(f'{file} {ID} {missing_count} {per_GC:.0f}')
        sys.stdout

        '''generate fasta file for the fastq file which will only contain reads/sequence with 
        quality score greater than the quality filter value provided on the command line'''
        
        if qscore > FastqParser.quality_filter: # if quality score is greater
            outfile = open(filepath, mode='a') # open a fasta file
            outfile.write('>' + ID) # append the ID of the sequence
            outfile.write('\n' + seq[1]) # append the sequence into the file