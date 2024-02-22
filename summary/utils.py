import re
import os
import pandas as pd
from math import log10
import matplotlib.pyplot as plt

def extractValueFromLine(lines, line, num, value_type):
    value = lines[lines.index(line)+num]
    if value_type == "int":
        return int(re.findall(r'\d+', value)[0])
    elif value_type == "float":
        # example: value = "thres = 1e-05" -> output 1e-05
        # 1. divide by = and remove ' '
        value = value.split('=')[1].strip()
        # 2. convert to float
        return float(value)
    elif value_type == "memory" or value_type == "memory_pick":
        memory, pick = getMemoryVals(value)
        if value_type == "memory":
            return memory
        else:
            return pick
    else:
        return None

def getResultSingleJob(job):
    output_all = getOutputFilenames(job)
    df1 = getResultStep1(output_all["step1"])
    df2 = getResultStep2(output_all["step2"])
    df3 = getResultStep3(output_all["step3"])
    df4 = getResultStep4(output_all["step4"])
    df5 = getResultStep5(output_all["step5"])
    df_total = pd.concat([df1, df2, df3, meanOfStep4(df4), df5], axis=1, join='outer')
    n_genes = df4['n_gene'].sum()
    n_batch = len(df4['4_num_batch'])
    df_param2 = pd.DataFrame({'job': [job], 'n_genes': [n_genes], 'num_batch': [n_batch]}, index=[0])
    df_total = pd.concat([df_total, df_param2], axis=1)
    return df_total

def getTotalResults(jobs_list):
    df_total = pd.DataFrame()
    for job in jobs_list:
        df = getResultSingleJob(job)
        df_total = pd.concat([df_total, df], axis=0)
    return df_total

def getOutputFilenames(job):
    OUTPUT_DIR = os.path.join(job, "slurm_out")
    OUTPUT_DIR_STEP4 = os.path.join(OUTPUT_DIR, "step4")
    OUTPUT_FILES_BATCH = os.listdir(OUTPUT_DIR_STEP4)
    output_file1 = os.path.join(OUTPUT_DIR, "step1.txt")
    output_file2 = os.path.join(OUTPUT_DIR, "step2.txt")
    output_file3 = os.path.join(OUTPUT_DIR, "step3.txt")
    output_file5 = os.path.join(OUTPUT_DIR, "step5.txt")
    output_file4 = [os.path.join(OUTPUT_DIR_STEP4, f) for f in OUTPUT_FILES_BATCH]
    output_all = {
        "step1": output_file1,
        "step2": output_file2,
        "step3": output_file3,
        "step4": output_file4,
        "step5": output_file5
    }
    return output_all

def getFileSizeFromLine(line):
    equal_sign_index = line.index('=')
    val_str = line[equal_sign_index + 1:].strip()
    if 'KB' in line:
        kb_index = val_str.index('KB')
        val = val_str[:kb_index].strip()
    else:
        gb_index = val_str.index('GB')
        val = val_str[:gb_index].strip()
    return val
    
def getMemoryVals(line):
    # Find the index of '(' and ')' to isolate the memory and pick memory information
    open_bracket_index = line.find('(')
    close_bracket_index = line.find(')')

    # Extract memory usage and pick memory substrings
    memory_info_pick = line[open_bracket_index + 1:close_bracket_index]

    # Split the memory info into memory and pick memory parts
    memory_info_final = line[:open_bracket_index]
    # Extract memory usage and pick memory values and remove leading and trailing spaces
    memory_usage = float(re.findall(r'\d+', memory_info_final)[0])
    pick_memory = float(re.findall(r'\d+', memory_info_pick)[0])
    
    # check if the memory is in GB or MB
    if 'MB' in memory_info_final:
        memory_usage = memory_usage / 1024.
    if 'MB' in memory_info_pick:
        pick_memory = pick_memory / 1024.
    return memory_usage, pick_memory


def getResultStep1(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    df = pd.DataFrame()
    row = {}
    for line in lines:
        if "default parameters" in line:
            row["k"] = extractValueFromLine(lines, line, 1, "int")
            row["thres"] = extractValueFromLine(lines, line, 2, "float")
            row["target"] = extractValueFromLine(lines, line, 3, "int")
        
        if " === Duration summaries ===" in line:
            row['1_readFastaFile'] = extractValueFromLine(lines, line, 1, "int")
            row['1_computeKmerTable'] = extractValueFromLine(lines, line, 2, "int")
            row['1_save'] = extractValueFromLine(lines, line, 3, "int")
            row['1_Memory usage'] = extractValueFromLine(lines, line, 7, "memory")
            row['1_pick Memory usage'] = extractValueFromLine(lines, line, 7, "memory_pick")
    df = df.append(row, ignore_index=True)
    return df

def getResultStep2(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    df = pd.DataFrame()
    row = {}
    for line in lines:
        if " === Duration summaries ===" in line:
            line_context = lines[lines.index(line) + 1]
            context = getFileSizeFromLine(lines[lines.index(line) + 5])
            public = getFileSizeFromLine(lines[lines.index(line) + 6])
            private = getFileSizeFromLine(lines[lines.index(line) + 7])
            mult = getFileSizeFromLine(lines[lines.index(line) + 8])
            row['2_GenCryptoContext'] = extractValueFromLine(lines, line, 1, "int")
            row['2_KeyGen'] = extractValueFromLine(lines, line, 2, "int")
            row['2_Serialize'] = extractValueFromLine(lines, line, 3, "int")
            row['2_context'] = context
            row['2_public'] = public
            row['2_private'] = private
            row['2_eval-mult'] = mult
            row['2_Memory usage'] = extractValueFromLine(lines, line, 9, "memory")
            row['2_Pick Memory usage'] = extractValueFromLine(lines, line, 9, "memory_pick")    
    df = df.append(row, ignore_index=True)
    return df



def getResultStep3(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    df = pd.DataFrame()
    row = {}
    for line in lines:
        if "n_kmer_total = " in line:
            row['n_kmer_total'] = extractValueFromLine(lines, line, 0, "int")
        if " === Duration summaries ===" in line:
            row['3_readFastaFile'] = extractValueFromLine(lines, line, 1, "int")
            row['3_load_kmerList'] = extractValueFromLine(lines, line, 2, "int")
            row['3_encodeRead'] = extractValueFromLine(lines, line, 3, "int")
            row['3_EncryptReadKmer'] = extractValueFromLine(lines, line, 4, "int")
            row['3_serializeCtxtRead'] = extractValueFromLine(lines, line, 5, "int")
            row['3_Memory usage'] = extractValueFromLine(lines, line, 11, "memory")
            row['3_Pick Memory usage'] = extractValueFromLine(lines, line, 11, "memory_pick")
    df = df.append(row, ignore_index=True)
    return df

def getResultStep4(files):
    df_total = pd.DataFrame()
    for filename in files:
        num_batch = int(re.findall(r'\d+', filename)[-1])
        with open(filename, 'r') as file:
            lines = file.readlines()
        row = {}

        for line in lines:
            if " == conditional parameters ==" in line:
                gene_start = extractValueFromLine(lines, line, 2, "int")
                gene_end = extractValueFromLine(lines, line, 3, "int")
                n_gene = gene_end - gene_start + 1
                row['n_gene'] = n_gene
                n_gene = float(n_gene)
            if " === Duration summaries ===" in line:
                row['4_loadKmerTable'] = extractValueFromLine(lines, line, 1, "int")
                row['4_loadKmerList'] = extractValueFromLine(lines, line, 2, "int")
                row['4_encodeReference'] = extractValueFromLine(lines, line, 3, "int") / n_gene
                row['4_readCtxtRead'] = extractValueFromLine(lines, line, 4, "int")
                row['4_multCtxtByEncodedRef'] = extractValueFromLine(lines, line, 5, "int") / n_gene
                row['4_serializeCtxtOut'] = extractValueFromLine(lines, line, 6, "int") / n_gene
                row['4_Memory usage'] = extractValueFromLine(lines, line, 11, "memory")
                row['4_Memory_per_gene'] = extractValueFromLine(lines, line, 11, "memory") / n_gene
                row['4_Pick Memory usage'] = extractValueFromLine(lines, line, 11, "memory_pick")
                row['4_num_batch'] = num_batch
        df_total = df_total.append(row, ignore_index=True)
    return df_total

def meanOfStep4(df4):
    # remove batch column
    df4 = df4.drop(columns=["4_num_batch"])
    # compute mean of each column and save as a one-row dataframe
    df4_mean = df4.mean()
    df4_mean = df4_mean.to_frame().transpose()
    return df4_mean

def getResultStep5(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    
    df = pd.DataFrame()
    row = {}
    for line in lines:
        if " === Duration summaries ===" in line:
            row['5_loadKmerTable'] = extractValueFromLine(lines, line, 1, "int")
            row['5_readCtxtOut'] = extractValueFromLine(lines, line, 2, "int")
            row['5_decryptOutput'] = extractValueFromLine(lines, line, 3, "int")
            row['5_Memory usage'] = extractValueFromLine(lines, line, 5, "memory")
            row['5_Pick Memory usage'] = extractValueFromLine(lines, line, 5, "memory_pick")
    df = df.append(row, ignore_index=True)
    return df