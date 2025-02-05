import re
import pandas as pd
from math import log10
import matplotlib.pyplot as plt
import os
import argparse

def extractValueFromLine(lines, line, num, value_type):
    line_num = lines.index(line) + num
    if line_num >= len(lines):
        line_num -= 1
    value = lines[line_num]
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

def getTotalResults(jobs_list):
    df_total = pd.DataFrame()
    for job in jobs_list:
        df = getResultSingleJob(job)
        df_total = pd.concat([df_total, df], axis=0)
    return df_total

def getOutputFilenames(job):
    OUTPUT_DIR = job + "/log/"
    OUTPUT_DIR_STEP4 = os.path.join(OUTPUT_DIR, "step4")
    OUTPUT_FILES_BATCH = os.listdir(OUTPUT_DIR_STEP4)
    output_file1 = os.path.join(OUTPUT_DIR, "step1_generate_kmerTable_cloud.txt")
    output_file2 = os.path.join(OUTPUT_DIR, "step2_he_keygen_local.txt")
    output_file3 = os.path.join(OUTPUT_DIR, "step3_encode_and_encrypt_local.txt")
    output_file5 = os.path.join(OUTPUT_DIR, "step5_decrypt_and_return_gene_vector_local.txt")
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


def getResultSingleJob(output_all):
    row = {}
    
    # Step 1: Read parameters and timings
    step1_file = output_all['step1']
    with open(step1_file, 'r') as file:
        lines = file.readlines()
    for line in lines:
        # if "default parameters" in line:
        #     row["k"] = extractValueFromLine(lines, line, 1, "int")
        #     row["thres"] = extractValueFromLine(lines, line, 2, "float")
        #     row["target"] = extractValueFromLine(lines, line, 3, "int")
        if " === Duration summaries ===" in line:
            row['1_readFastaFile'] = extractValueFromLine(lines, line, 1, "int")
            row['1_computeKmerTable'] = extractValueFromLine(lines, line, 2, "int")
            row['1_save'] = extractValueFromLine(lines, line, 3, "int")
            row['1_memory'] = extractValueFromLine(lines, line, 7, "memory")
            row['1_memory_pick'] = extractValueFromLine(lines, line, 7, "memory_pick")
        
    
    # Step 2: Read timings
    step2_file = output_all['step2']
    with open(step2_file, 'r') as file:
        lines = file.readlines()
    for line in lines:
        if " === Duration summaries ===" in line:
            row['2_GenCryptoContext'] = extractValueFromLine(lines, line, 1, "int")
            row['2_KeyGen'] = extractValueFromLine(lines, line, 2, "int")
            row['2_Serialize'] = extractValueFromLine(lines, line, 3, "int")
            row['2_memory'] = extractValueFromLine(lines, line, 9, "memory")
            row['2_memory_pick'] = extractValueFromLine(lines, line, 9, "memory_pick")
    
    # Step 3: Read timings
    step3_file = output_all['step3']
    with open(step3_file, 'r') as file:
        lines = file.readlines()
    for line in lines:
        if " === Duration summaries ===" in line:
            row['3_readFastaFile'] = extractValueFromLine(lines, line, 1, "int")
            row['3_load_kmerList'] = extractValueFromLine(lines, line, 2, "int")
            row['3_encodeRead'] = extractValueFromLine(lines, line, 3, "int")
            row['3_EncryptReadKmer'] = extractValueFromLine(lines, line, 4, "int")
            row['3_serializeCtxtRead'] = extractValueFromLine(lines, line, 5, "int")
            row['3_memory'] = extractValueFromLine(lines, line, 11, "memory")
            row['3_memory_pick'] = extractValueFromLine(lines, line, 11, "memory_pick")
    
    # Step 4: Read timings for all batches
    step4_files = output_all['step4']
    # Initialize lists to store values for each batch
    loadKmerTable_times = []
    loadKmerList_times = []
    readCtxtRead_times = []
    encodeAndMult_times = []
    serializeCtxtOut_times = []
    memory_usage = []
    memory_pick = []
    n_gene_total = 0
    encodeAndMult_time_max = 0
    encodeAndMult_time_min = 2**63-1

    for filename in step4_files:
        with open(filename, 'r') as file:
            lines = file.readlines()
        for line in lines:
            if "n_genes = " in line:
                n_gene_total += extractValueFromLine(lines, line, 0, "int")
            if " === Duration summaries ===" in line:
                loadKmerTable_times.append(extractValueFromLine(lines, line, 1, "int"))
                loadKmerList_times.append(extractValueFromLine(lines, line, 2, "int"))
                readCtxtRead_times.append(extractValueFromLine(lines, line, 3, "int"))
                encodeAndMult_current = extractValueFromLine(lines, line, 4, "int")
                encodeAndMult_times.append(encodeAndMult_current)
                if encodeAndMult_current > encodeAndMult_time_max:
                    encodeAndMult_time_max = encodeAndMult_current
                if encodeAndMult_current < encodeAndMult_time_min:
                    encodeAndMult_time_min = encodeAndMult_current
                
                serializeCtxtOut_times.append(extractValueFromLine(lines, line, 5, "int"))
                memory_usage.append(extractValueFromLine(lines, line, 10, "memory"))
                memory_pick.append(extractValueFromLine(lines, line, 10, "memory_pick"))
    
    print("n_gene_total", n_gene_total)
    # Calculate mean values and add to the row
    batch_num = len(loadKmerTable_times)
    assert batch_num == len(loadKmerList_times)
    assert batch_num == len(readCtxtRead_times)
    assert batch_num == len(encodeAndMult_times)
    assert batch_num == len(serializeCtxtOut_times)
    # row['4_loadKmerTable_mean'] = sum(loadKmerTable_times) / len(loadKmerTable_times)
    # row['4_loadKmerList_mean'] = sum(loadKmerList_times) / len(loadKmerList_times)
    # row['4_readCtxtRead_mean'] = sum(readCtxtRead_times) / len(readCtxtRead_times)
    # row['4_encodeAndMult_mean'] = sum(encodeAndMult_times) / len(encodeAndMult_times)
    # row['4_serializeCtxtOut_mean'] = sum(serializeCtxtOut_times) / len(serializeCtxtOut_times)
    # row['4_memory'] = sum(memory_usage) / len(memory_usage)
    # row['4_memory_pick'] = sum(memory_pick) / len(memory_pick)
    row['4_loadKmerTable_mean'] = sum(loadKmerTable_times) / batch_num
    row['4_loadKmerList_mean'] = sum(loadKmerList_times) / batch_num
    row['4_readCtxtRead_mean'] = sum(readCtxtRead_times) / batch_num
    row['4_encodeAndMult_mean'] = sum(encodeAndMult_times) / batch_num
    row['4_encodeAndMult_max'] = max(encodeAndMult_times)
    row['4_encodeAndMult_min'] = min(encodeAndMult_times)
    row['4_encodeAndMult'] = encodeAndMult_times
    row['4_encodeAndMult_std'] = (sum([(x - sum(encodeAndMult_times) / batch_num) ** 2 for x in encodeAndMult_times]) / batch_num) ** 0.5
    row['4_serializeCtxtOut_mean'] = sum(serializeCtxtOut_times) / batch_num
    row['4_memory'] = sum(memory_usage) / batch_num
    row['4_memory_pick'] = sum(memory_pick) / batch_num
    
    # Step 5: Read timings
    step5_file = output_all['step5']
    with open(step5_file, 'r') as file:
        lines = file.readlines()
    for line in lines:
        if " === Duration summaries ===" in line:
            row['5_loadKmerTable'] = extractValueFromLine(lines, line, 1, "int")
            row['5_readCtxtOut'] = extractValueFromLine(lines, line, 2, "int")
            row['5_decryptOutput'] = extractValueFromLine(lines, line, 3, "int")
            row['5_memory'] = extractValueFromLine(lines, line, 5, "memory")
            row['5_memory_pick'] = extractValueFromLine(lines, line, 5, "memory_pick")
    
    # Convert the row dictionary to a DataFrame with a single row
    df = pd.DataFrame([row])
    return df


# BASE_FOLDER = "out/"
# job = BASE_FOLDER + "20240927_140838"
# print("reading ", job)
# outputs = getOutputFilenames(job)
# df = getResultSingleJob(outputs)
# print(df)
# df = df.transpose()
# df.columns = ["time (ms)"]
# df.to_csv(job + "/time_summary.csv", index=True, index_label="job")

if __name__ == "__main__":
    print("run")
    parser = argparse.ArgumentParser(description='Process job outputs and summarize timings.')
    parser.add_argument('job_id', type=str, help='The job ID to process')
    args = parser.parse_args()
    job_id = args.job_id
    BASE_FOLDER = ""
    location = BASE_FOLDER + job_id
    print("dir ", location)
    outputs = getOutputFilenames(location)
    df = getResultSingleJob(outputs)
    
    print(df)
    df = df.transpose()
    df.columns = ["time (ms)"]
    df.to_csv(location + "/time_summary.csv", index=True, index_label="job")
