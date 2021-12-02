import io
import re
import sys
import time

'''Part 1'''

'''A function to count kmers using lists'''

def count_kmers_non_dict(seqs, k):
    
    # seqs should be input as a list
    
    # initializing the lists
    kmer_ls = list()
    unique_kmer = list()
    
    # default valid characters
    valid = ["A", "T", "C", "G", "a", "t", "c", "g"]
    
    # for each strings in the list of strings
    for seq in seqs:
        # iterating through the sequence:
        for i in range(len(seq)-k+1):

            # obtain the kmer for each iteration
            kmer = seq[i:i+k]

            # see if the kmer contains invalid characters
            valid_kmer = True
            for char in kmer:
                if char not in valid:
                    valid_kmer = False
                    break

            kmer = kmer.upper()

            # appending the list for valid kmers
            if valid_kmer == True:
                kmer_ls.append(kmer)

                # to create a unique list
                if kmer not in unique_kmer:
                    unique_kmer.append(kmer)

    # Sorting the list for a better output presentation
    unique_kmer = sorted(unique_kmer)
    
    # output
    kmer_count = []
    for kmer in unique_kmer:
        kmer_count.append(kmer_ls.count(kmer))
    
    return unique_kmer, kmer_count


'''A function to return the output'''

def write2file_non_dict(kmer_list, counts, outfilename):
    
    # open the output file
    with open(outfilename, "w") as fh:
        
        # write the kmer and its count line by line
        for i in range(len(kmer_list)):
            fh.write(f"{kmer_list[i]} : {counts[i]}\n")
            
            
'''Part 2'''

'''A function to count kmers using dicts'''

def count_kmers(seqs, k):
    
    # initializing the dict
    dc = dict()
    
    # default valid characters
    valid = ["A", "T", "C", "G", "a", "t", "c", "g"]
    
    # for each string in the list of strings:
    for seq in seqs:
        
        # iterating through the sequence:
        for i in range(len(seq)-k+1):

            # obtain the kmer for each iteration
            kmer = seq[i:i+k]

            # see if the kmer contains invalid characters
            valid_kmer = True
            for char in kmer:
                if char not in valid:
                    valid_kmer = False
                    break

            kmer = kmer.upper()

            # appending the list for valid kmers
            if valid_kmer == True:
                if kmer in dc.keys():
                    dc[kmer] += 1
                else:
                    dc[kmer] = 1
                
    # Output
    return dc

'''A function to write the output of the kmer dict'''

def write2file(kmer_dict, outfilename):
    
    # open the file
    with open(outfilename, "w") as fh:
    
        # writing line by line for each key in dict and its value
        for key in sorted(kmer_dict.keys()):
            fh.write(f"{key} : {kmer_dict[key]}\n")
            
"""A function to read .txt sequence file"""

def r_txt(filepath):
    
    # number of lines in the file
    no = sum(1 for i in io.open(filepath))
    
    # initialization
    out = ""
    
    with io.open(filepath) as fh:
        
        for i in range(no):
            line = fh.readline().strip()
            out += line
            
    return [out]


'''Part 3'''

"""A function to read fasta / fa sequence file"""

def r_fasta(filepath):
    
    # number of lines in the file
    no = sum(1 for i in io.open(filepath))
    
    # initialization
    out = []
    seq = ""
    new = True
    
    with io.open(filepath) as fh:
        
        for i in range(no):
            
            line = fh.readline().strip()
            
            # for each header
            if re.search(r">", line):
                new = True
                if len(seq) != 0:
                    out.append(seq)
                
            # for the sequence lines
            else:
                
                # second line onwards
                if new == False:
                    seq += line
                    # for the very last line in the file
                    if i == (no-1):
                        out.append(seq)
                        
                # first line of the sequence
                else:
                    
                    new = False
                    seq = line
            
    return out


"""A function to read fastq sequence file"""

def r_fastq(filepath):
    
    # number of lines in the file
    no = sum(1 for i in io.open(filepath))
    
    # initialization
    out = []
    seq = ""
    seq_line = 0
    
    with io.open(filepath) as fh:
        
        for i in range(no):
            
            line = fh.readline().strip()
            
            # The first line in the 4-line fastq format; locating the second line (real sequence) as seq_line
            if re.search(r"@", line):
                
                # getting the last digit of the label line
                digit = line[-1]
                
                # this is for indicating the computer to add the sequence for the next iteration
                seq_line = i+1
                
            # Adding the second line in the 4-line fastq format
            if seq_line == i:
                
                if (digit == "1") and (len(seq) != 0):
                    out.append(seq)
                    seq = line
                else:
                    seq += line
                
            if i == (no-1):
                out.append(seq)
            
    return out



'''Main program'''

def main():
    
    
    ### Part 1
    
    # calculate and output the kmer counts using the list method

    seq = ['ctccaaagaaattgtagttttcttctggcttagaggtagatcatcttggtccaatcagactgaaatgccttgaggctagatttcagtctttgtGGCAGCTGgtgaatttctagtttgccttttcagctagggattagctttttaggggtcccaatgcctagggagatttctaggtcctctgttccttgctgacctccaat']
    kmer, count = count_kmers_non_dict(seq, 2)
    write2file_non_dict(kmer, count, "part1_out.txt")
    
    
    ### Part 2
    
    # calculate and output the kmer counts using the dict method
    kmer_dc = count_kmers(seq,2)
    write2file(kmer_dc, "part2_out.txt")
    
    # (Part 1 and part 2 will be generated automatically)
    
    
    ### Part 3
    
    # input from the command line
    filename = sys.argv[1]
    extension = filename.split(".")[-1]
    ks = [int(i) for i in sys.argv[2].split(",")]
    try:
        if sys.argv[3] == "-t":
            option = sys.argv[3]
    except:
        option = "-d"
    
    # decision made from sys.argv[3] (-t)
    if option == "-t":

        # comparing the time used for method 1 and 2
        seq = r_txt(filename)
        with open("readme.txt", "w") as fh:

            # First method
            start = time.time()
            count_kmers_non_dict(seq, ks[0])
            end = time.time()
            fh.write(f"Time elapsed for the first method is {end-start} seconds.\n")

            # Second method
            start = time.time()
            count_kmers(seq, ks[0])
            end = time.time()
            fh.write(f"Time elapsed for the second method is {end-start} seconds.\n")
            
    else:
        
        # recognizing the file format and extract the sequence
        if (extension == "fasta") | (extension == "fa") | (extension == "faa") | (extension == "fna"):
            seq = r_fasta(filename)
        elif (extension == "fastq") | (extension == "fq"):
            seq = r_fastq(filename)
        elif (extension == "txt"):
            seq = r_txt(filename)
        else:
            print("file extention not supported")
            sys.exit()
        
        # for each kmer in the kmer list input in the command line
        for kk in ks:
            kmer_dc = count_kmers(seq,kk)
            write2file(kmer_dc, "part3_out_" + str(kk) + ".txt")
            

# Running the program

main()
    