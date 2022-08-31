#!/usr/bin/env python

#Bash commands used for Part 1 to process STAR alignment file. See Alignment.srun and Genome_gen.srun for STAR slurm scripts.
#See labbook_PS8.py for full command line process.
    #samtools view -S -b STAR_align_dre_WT_ovar12Aligned.out.sam > STAR_align_dre_WT_ovar12Aligned.out.bam
    #samtools view STAR_align_dre_WT_ovar12Aligned.out.bam | head
    #samtools sort STAR_align_dre_WT_ovar12Aligned.out.bam -o STAR_align_dre_WT_ovar12Aligned.sorted.bam
    #samtools index STAR_align_dre_WT_ovar12Aligned.sorted.bam
    #samtools view STAR_align_dre_WT_ovar12Aligned.sorted.bam 1 > STAR_align_dre_WT_ovar12Aligned_chr1.sam
    #wc -l STAR_align_dre_WT_ovar12Aligned_chr1.sa #To count number of records in chr1 file

#file = "STAR_align_dre_WT_ovar12Aligned.out.sam" - PS8
#file = "test.sam" - PS8
#output = "PS8_output.txt" - PS8

#Data files QAA
file = "/projects/bgmp/shared/2017_sequencing/demultiplexed/14_3B_control_S10_L008_R1_001.fastq.gz"
file = "/projects/bgmp/shared/2017_sequencing/demultiplexed/14_3B_control_S10_L008_R2_001.fastq.gz"
file = "/projects/bgmp/shared/2017_sequencing/demultiplexed/11_2H_both_S9_L008_R1_001.fastq.gz"
file = "/projects/bgmp/shared/2017_sequencing/demultiplexed/11_2H_both_S9_L008_R2_001.fastq.gz"

#Output file QAA
output = "QAA_output.txt"

#Initialize lists for mapped reads and unmapped reads, and variable for line/record counter
reads_dict = {}

#Open SAM input file from STAR alignment
with open(file,"r") as fh:
    while True:
        #Strip newlines
        record = fh.readline().strip()
        if record == "":
            #EOF
            break
        if record.startswith("@") == False:
            #For each record (line) in file that is not a header line
            flag = int(record.split("\t")[1]) #Assign bitwise flag as int to flag variable
            name = record.split("\t")[0] #Name of record
            reads_dict[name] = flag #Create dictionary object (name=key,flag=value)

    #Dictionaries have unique keys, so multiple alignment reads will be ignored

    #Initialize counters for mapped and unmapped reads
    mapped_count = 0
    unmapped_count = 0

    #For each read (key) in dictionary
    for read in reads_dict:
        #Check bitwise flag (value) for if read is mapped/primary alignment:
        if((reads_dict[read] & 4) != 4):
            if ((reads_dict[read] & 256 ) != 256):
                #If mapped and primary alignment, increase mapped_count
                mapped_count += 1
        else:
            #Else, increase unmapped_count
            unmapped_count += 1
    
    #Sum mapped and unmapped counts for total count
    total_count = mapped_count + unmapped_count

    #Print all three counts (mapped, unmapped, total count)
    print("Mapped Count: " + str(mapped_count))
    print("Unmapped Count: " + str(unmapped_count))
    print("Total Read Count: " + str(total_count))

    #Write counts to output (txt) file
    with open(output,"w") as fh2:
        fh2.write("Mapped Count: "+str(mapped_count)+"\n"+"Unmapped Count: "+str(unmapped_count)+"\n"+"Total Read Count: "+str(total_count))
    
    #Mapped Count: 
    #Unmapped Count: 
    #Total Read Count: 
