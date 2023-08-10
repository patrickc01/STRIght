import os
import glob
import sys
import csv
import numpy as np
from collections import Counter
import pandas as pd 

#This script analyzes an NGS run for STRs.  STRs are loaded from the file STR_data_v14.csv.
#Next the script goes through each fastq folder in the targe directory.  Within the fastq file, each line is checked for being a STR.  If an STR is found, it is recorded 
#as an STR within that fastq file.  After iterating through all lines of the fastq file, the top 4 reads are found and recorded.  Then the script moves on to the next fastq file.

class well():
    #Well represents a fastq file containing NGS reads for a specific sample ie human_sample_3.
    #This fastq file will be scanned for all possible STR present in the STR data sheet.
    def __init__(self,fastq_file):
        self.well_name = fastq_file
        self.hits = []          #list of str_object hits found in th well
        self.STR_library = []   #list of STR objects for the well
        self.load_str_objects()     
        self.plate_results=pd.DataFrame() 
        self.all_hits_list=[]

    def append_hit(self, STR_hit):
        self.hits.append(STR_hit)

    def load_str_objects(self):
        with open("STRight_key.csv", 'r') as csvfile:
            str_reader = csv.reader(csvfile)
            next(str_reader, None)
            for row in str_reader:
                self.STR_library.append(make_str_objects(row))

    def condense_hits(self):# STR_library):
        #Goes through each fastq file and searches for all STRs found in the fastq file, condenses results to top 4 reads     
        #Data used to generate the .csv file   
        top_number_of_hits=4
        for x in self.STR_library:
            if x.c_counter > 100:          #There must be at least 90 reads for the STR to be considered valid
                total_count = x.c_counter
                x.top_reads = Counter(x.STR_lengths_list).most_common(top_number_of_hits)
                x.SNP_test = format(x.start_counter / x.end_counter, '.2f')
                unique_STR_list=[self.well_name,x.STR_name, x.c_counter, x.SNP_test]
                for top_STR_length,top_STR_count in x.top_reads:
                    try:
                        ratio= top_STR_count/total_count *100
                        unique_STR_list.extend([top_STR_length,format(ratio, '.2f')])
                    except IndexError:
                        pass
                self.all_hits_list.append(unique_STR_list)
            else:
                pass
        return self.all_hits_list

class make_str_objects():
    def __init__(self,row):
        self.STR_name = row[0]
        self.sequence = str.upper(row[1])
        self.start = str.upper(row[2])
        self.end = str.upper(row[3])
        self.STR_size= int(row[4])          #Size of STR, usually 4bp or 5bp        
        self.BP_modifier= int(row[5])         #modifier, used to make control samples match known STR repeat lengths in complex cases.  modifier is an integer of repeat #
        self.ref_repeat_num=row[6]         # repeats in reference STR_data.csv file
        self.SNP_mod = str(row[7])
        self.normal_distance = self.sequence.find(self.end) - self.sequence.find(self.start)-len(self.start)
        self.SNP_test = []
        self.c_counter = 0
        self.start_counter=0        #Number of times start sequence is found in a fastq file
        self.end_counter = 0        #Number of times end sequences if found in a fastq file
        self.STR_lengths_list = []  #List of STR_lengths found in each line of fastq file
        self.STR_sequences_list= [] #List of all sequences found for each unique STR in a fastq file, used to find most common sequences
        self.top_reads = []         #top 2 read counts for STR
        self.test_mod_repeats=(self.normal_distance-self.BP_modifier) /self.STR_size
        self.condensed_list=[]
        try: 
            self.SNP_mod_dict = dict(map(lambda x: x.split('='), self.SNP_mod.split(', ')))
            for key, value in self.SNP_mod_dict.items():
                self.SNP_mod_dict[key]=int(value)
        except ValueError:
            self.SNP_mod_dict=False

    def append_length(self,STR_size):
        self.STR_lengths_list.append(STR_size)

def write_to_file(record_entry,file_name):
    #Write STR sequences and counts to .txt file
    f = open(file_name, "w")
    for each_well in record_entry:
        for element in each_well:
            print(element, file=f)
        f.write('\n')
    f.close()

def analyze_well(well_name):
    well_fastq = open(str(well_name), "r")     #open the current fastq file
    well_class = well(well_name)               #Make a well class that contains all potential STRs, with the name being the name of the well
    for line in well_fastq:                    #Go through each line of the fastq
        for individual_STR in well_class.STR_library:                                   #test each line for every STR in STR_libary
            if line.find(individual_STR.start)>0 and line.find(individual_STR.end)>0:   #If line looks like an STR, record STR size
                start = line.find(individual_STR.start)
                end = line.find(individual_STR.end)+len(individual_STR.end)
                individual_STR.start_counter +=1
                individual_STR.end_counter +=1
                individual_STR.c_counter += 1
                #Check if STR specific SNPs are found in sequence and modify length if needed
                STR_size = (line.find(individual_STR.end)- line.find(individual_STR.start)-len(individual_STR.start))- individual_STR.BP_modifier
                individual_STR.STR_sequences_list.append(line[start:end+9]+" ("+ str(STR_size/individual_STR.STR_size)+")") #append the sequence of the found STR repeat
                curr_length = (STR_size)//individual_STR.STR_size + (STR_size % individual_STR.STR_size/10)
                
                if individual_STR.SNP_mod_dict:
                    for key in individual_STR.SNP_mod_dict:
                        if key in line:
                            curr_length += individual_STR.SNP_mod_dict[key] #append the size of the found STR repeat
                individual_STR.append_length(curr_length)
            elif line.find(individual_STR.start)>0 and line.find(individual_STR.end)<0:  #If STR start sequence is found, but end sequence is not found, only increment start counter.  Used for SNP test
                individual_STR.start_counter +=1
            elif line.find(individual_STR.end)>0 and line.find(individual_STR.start)<0:  #If STR start sequence is found, but end sequence is not found, only increment start counter.  Used for SNP test
                individual_STR.end_counter +=1
            else:
                pass
    return well_class

def search_NGS_run():
    NGS_all_wells = []        #Master list of all samples from the NGS run, each sample consists of a well Class, each well class contains potential STRs
    fastq_files = '*.fastq'   #fastq files to analyze
    file_name = 'STR_analysis.txt'    #Name of output txt file
    master_Record = []
    
    print("Program Running")
    for well_name in glob.glob(fastq_files):            #Analyze every well (fastq file) for all possible STRs
        NGS_all_wells.append(analyze_well(well_name))   #Create a list of all wells, each well will contain information on all STRs
    for each_well in NGS_all_wells:
        print('checking well {} for hits'.format(each_well.well_name))
        for individual_STR in each_well.STR_library:
            if individual_STR.c_counter >100:
                top_4_reads = Counter(individual_STR.STR_sequences_list).most_common(4)
                record1 = ([each_well.well_name + " " + "TOTAL:" + str(individual_STR.c_counter) + " " + "STR:" + individual_STR.STR_name])
                for k, v in top_4_reads:
                    record1.append(str('{}, {}'.format(k,v)))   #Append the top 6 reads found in the NGS file to the analysis file
                master_Record.append(record1)
            else:
                pass
    master_Record = sorted(master_Record)
    write_to_file(master_Record,file_name)  #Write .txt file containing top read sequences for each well/STR
    final_csv_list = []
    for each_well in NGS_all_wells:
        condensed_STR_list = each_well.condense_hits()#STR_library)
        if condensed_STR_list:
            final_csv_list.extend(condensed_STR_list)
        else:
            pass
    final_df=pd.DataFrame(final_csv_list, columns=['Well','STR','Total_count','SNP_test','STR_length_1','Fraction_1','STR_length_2','Fraction_2','STR_length_3','Fraction_3','STR_length_4','Fraction_4'])
    final_df=final_df.sort_values(by=['Well','STR'])
    final_df.to_csv('STR_results_df.csv', index=False)

if __name__=="__main__":
    search_NGS_run()
    print("Done")
