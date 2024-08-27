import pandas as pd
import os
import shutil
import pyhmmer
import subprocess
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
warnings.filterwarnings("ignore", category=UserWarning, module="scipy")
from Bio import SeqIO
from pyhmmer.easel import TextSequence, TextMSA, DigitalSequenceBlock
import warnings
warnings.filterwarnings("ignore", category=UserWarning, module="numpy")
warnings.filterwarnings("ignore", category=UserWarning, module="scipy")
import glob
import csv
import matplotlib.pyplot as plt

# L2 名稱集合
L2_name = {
    "sp|Q9NQ29|LUC7L_HUMAN", # 1
    "sp|Q9CYI4|LUC7L_MOUSE", # 2
    "tr|A0A3Q2UH51|A0A3Q2UH51_CHICK", # 3
    "tr|Q6GLC1|Q6GLC1_XENTR", # 4
    "tr|H3AZN1|H3AZN1_LATCH", # 5
    "tr|Q6IQD4|Q6IQD4_DANRE", # 6
    "sp|Q9Y383|LC7L2_HUMAN", # 7
    "sp|Q7TNC4|LC7L2_MOUSE", # 8
    "tr|Q5ZLW7|Q5ZLW7_CHICK", # 9
    "tr|Q28EN5|Q28EN5_XENTR", # 10
    "tr|H2ZXS0|H2ZXS0_LATCH", # 11
    "tr|A0A8M1NDA7|A0A8M1NDA7_DANRE", # 12
    "tr|A0A8C4NH22|A0A8C4NH22_EPTBU", # 13
    "tr|F6S0J6|F6S0J6_CIOIN", # 14
    "tr|A0A182H045|A0A182H045_AEDAL", # 15
    "tr|Q9VVI1|Q9VVI1_DROME", # 16
    "tr|B3RNA2|B3RNA2_TRIAD", # 17
    "tr|T2M333|T2M333_HYDVU", # 18
    "sp|Q09217|YP68_CAEEL", # 19
    "tr|A0A1X7V7G6|A0A1X7V7G6_AMPQE", # 20

    "tr|F4IZU3|F4IZU3_ARATH", # 21
    "tr|Q940U9|Q940U9_ARATH", # 22
    "tr|A0A3Q7GLL9|A0A3Q7GLL9_SOLLC", # 23
    "tr|U5D3Z2|U5D3Z2_AMBTC", # 24
    "tr|Q75LD6|Q75LD6_ORYSJ", # 25
    "tr|B4FUS0|B4FUS0_MAIZE", # 26
    "tr|B4FXT5|B4FXT5_MAIZE", # 27
    "tr|Q6K8R4|Q6K8R4_ORYSJ", # 28
    "tr|A0A2R6WMH7|A0A2R6WMH7_MARPO", # 29
    "tr|A0A388LKU4|A0A388LKU4_CHABU", # 30

    "tr|Q7S615|Q7S615_NEUCR", # 31
    "tr|Q5BCF4|Q5BCF4_EMENI", # 32
    "tr|F9X0V8|F9X0V8_ZYMTI", # 33
    "sp|Q07508|LUC7_YEAST", # 34
    "sp|Q9USM4|LUC7_SCHPO", # 35
    "tr|Q5KBY8|Q5KBY8_CRYNJ", # 36
    "tr|E3L016|E3L016_PUCGT", # 37
    "tr|A0A139ATP5|A0A139ATP5_GONPJ" # 38
}


# L3 名稱集合
L3_name = {
    "sp|Q5SUF2|LC7L3_MOUSE", # 1
    "sp|O95232|LC7L3_HUMAN", # 2
    "tr|F1NXR8|F1NXR8_CHICK", # 3
    "tr|Q28G85|Q28G85_XENTR", # 4
    "tr|H3ABC3|H3ABC3_LATCH", # 5
    "tr|Q1ED13|Q1ED13_DANRE", # 6
    "tr|B3RQ49|B3RQ49_TRIAD", # 7
    "tr|A0A8B6XHB2|A0A8B6XHB2_HYDVU", # 8
    "tr|F6Z8S6|F6Z8S6_CIOIN", # 9
    "tr|A0A8C4N567|A0A8C4N567_EPTBU", # 10
    "tr|A0A182GZD3|A0A182GZD3_AEDAL", # 11
    "tr|Q9W3X8|Q9W3X8_DROME", # 12
    "tr|Q8ITY5|Q8ITY5_CAEEL", # 13
    "tr|H2L0S8|H2L0S8_CAEEL", # 14
    "tr|A0A1X7UMY2|A0A1X7UMY2_AMPQE", # 15

    "tr|P94088|P94088_ARATH", # 16
    "tr|A0A3Q7GAS3|A0A3Q7GAS3_SOLLC", # 17
    "tr|W1P1B9|W1P1B9_AMBTC", # 18
    "tr|B4FE46|B4FE46_MAIZE", # 19
    "tr|Q84YS0|Q84YS0_ORYSJ", # 20
    "tr|A0A176W202|A0A176W202_MARPO", # 21
    "tr|A0A388KYX9|A0A388KYX9_CHABU", # 22
    "tr|A0A2R6XJM8|A0A2R6XJM8_MARPO", # 23
    "tr|A0A2K3CY04|A0A2K3CY04_CHLRE", # 24
    "tr|A8HP50|A8HP50_CHLRE", # 25
    "tr|A4SBD0|A4SBD0_OSTLU" # 26
}


classification = {
    "sp|Q9NQ29|LUC7L_HUMAN": 'Animal', # 1
    "sp|Q9CYI4|LUC7L_MOUSE": 'Animal', # 2
    "tr|A0A3Q2UH51|A0A3Q2UH51_CHICK": 'Animal', # 3
    "tr|Q6GLC1|Q6GLC1_XENTR": 'Animal', # 4
    "tr|H3AZN1|H3AZN1_LATCH": 'Animal', # 5
    "tr|Q6IQD4|Q6IQD4_DANRE": 'Animal', # 6
    "sp|Q9Y383|LC7L2_HUMAN": 'Animal', # 7
    "sp|Q7TNC4|LC7L2_MOUSE": 'Animal', # 8
    "tr|Q5ZLW7|Q5ZLW7_CHICK": 'Animal', # 9
    "tr|Q28EN5|Q28EN5_XENTR": 'Animal', # 10
    "tr|H2ZXS0|H2ZXS0_LATCH": 'Animal', # 11
    "tr|A0A8M1NDA7|A0A8M1NDA7_DANRE": 'Animal', # 12
    "tr|A0A8C4NH22|A0A8C4NH22_EPTBU": 'Animal', # 13
    "tr|F6S0J6|F6S0J6_CIOIN": 'Animal', # 14
    "tr|A0A182H045|A0A182H045_AEDAL": 'Animal', # 15
    "tr|Q9VVI1|Q9VVI1_DROME": 'Animal', # 16
    "tr|B3RNA2|B3RNA2_TRIAD": 'Animal', # 17
    "tr|T2M333|T2M333_HYDVU": 'Animal', # 18
    "sp|Q09217|YP68_CAEEL": 'Animal', # 19
    "tr|A0A1X7V7G6|A0A1X7V7G6_AMPQE": 'Animal', # 20

    "tr|F4IZU3|F4IZU3_ARATH": 'Plant', # 21
    "tr|Q940U9|Q940U9_ARATH": 'Plant', # 22
    "tr|A0A3Q7GLL9|A0A3Q7GLL9_SOLLC": 'Plant', # 23
    "tr|U5D3Z2|U5D3Z2_AMBTC": 'Plant', # 24
    "tr|Q75LD6|Q75LD6_ORYSJ": 'Plant', # 25
    "tr|B4FUS0|B4FUS0_MAIZE": 'Plant', # 26
    "tr|B4FXT5|B4FXT5_MAIZE": 'Plant', # 27
    "tr|Q6K8R4|Q6K8R4_ORYSJ": 'Plant', # 28
    "tr|A0A2R6WMH7|A0A2R6WMH7_MARPO": 'Plant', # 29
    "tr|A0A388LKU4|A0A388LKU4_CHABU": 'Plant', # 30

    "tr|Q7S615|Q7S615_NEUCR": 'Fungi', # 31
    "tr|Q5BCF4|Q5BCF4_EMENI": 'Fungi', # 32
    "tr|F9X0V8|F9X0V8_ZYMTI": 'Fungi', # 33
    "sp|Q07508|LUC7_YEAST": 'Fungi', # 34
    "sp|Q9USM4|LUC7_SCHPO": 'Fungi', # 35
    "tr|Q5KBY8|Q5KBY8_CRYNJ": 'Fungi', # 36
    "tr|E3L016|E3L016_PUCGT": 'Fungi', # 37
    "tr|A0A139ATP5|A0A139ATP5_GONPJ": 'Fungi', # 38

    #################################################

    "sp|Q5SUF2|LC7L3_MOUSE": 'Animal', # 1
    "sp|O95232|LC7L3_HUMAN": 'Animal', # 2
    "tr|F1NXR8|F1NXR8_CHICK": 'Animal', # 3
    "tr|Q28G85|Q28G85_XENTR": 'Animal', # 4
    "tr|H3ABC3|H3ABC3_LATCH": 'Animal', # 5
    "tr|Q1ED13|Q1ED13_DANRE": 'Animal', # 6
    "tr|B3RQ49|B3RQ49_TRIAD": 'Animal', # 7
    "tr|A0A8B6XHB2|A0A8B6XHB2_HYDVU": 'Animal', # 8
    "tr|F6Z8S6|F6Z8S6_CIOIN": 'Animal', # 9
    "tr|A0A8C4N567|A0A8C4N567_EPTBU": 'Animal', # 10
    "tr|A0A182GZD3|A0A182GZD3_AEDAL": 'Animal', # 11
    "tr|Q9W3X8|Q9W3X8_DROME": 'Animal', # 12
    "tr|Q8ITY5|Q8ITY5_CAEEL": 'Animal', # 13
    "tr|H2L0S8|H2L0S8_CAEEL": 'Animal', # 14
    "tr|A0A1X7UMY2|A0A1X7UMY2_AMPQE": 'Animal', # 15

    "tr|P94088|P94088_ARATH": 'Plant', # 16
    "tr|A0A3Q7GAS3|A0A3Q7GAS3_SOLLC": 'Plant', # 17
    "tr|W1P1B9|W1P1B9_AMBTC": 'Plant', # 18
    "tr|B4FE46|B4FE46_MAIZE": 'Plant', # 19
    "tr|Q84YS0|Q84YS0_ORYSJ": 'Plant', # 20
    "tr|A0A176W202|A0A176W202_MARPO": 'Plant', # 21
    "tr|A0A388KYX9|A0A388KYX9_CHABU": 'Plant', # 22
    "tr|A0A2R6XJM8|A0A2R6XJM8_MARPO": 'Plant', # 23
    "tr|A0A2K3CY04|A0A2K3CY04_CHLRE": 'Plant', # 24
    "tr|A8HP50|A8HP50_CHLRE": 'Plant', # 25
    "tr|A4SBD0|A4SBD0_OSTLU": 'Plant' # 26

}


def setup_directory(directory):
    """Ensure the directory exists."""
    if not os.path.exists(directory):
        os.makedirs(directory)

def classify_proteins(df, thresholds):
    """Classify proteins based on SPS thresholds."""
    df['Classification'] = df['SPS'].apply(
        lambda x: 'L2-type' if x > thresholds['L2'] else ('L3-type' if x < thresholds['L3'] else 'Unclassified')
    )
    return df

def run_external_command(cmd):
    """Run external shell commands."""
    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Failed to execute command: {cmd}")
        print(f"Error: {e}")

def calculate_sps(df):
    """Calculate the Subfamily Probability Score (SPS)."""
    df['SPS'] = (pd.to_numeric(df['L2'], errors='coerce') / pd.to_numeric(df['Length'], errors='coerce')) - (pd.to_numeric(df['L3'], errors='coerce') / pd.to_numeric(df['Length'], errors='coerce'))
    return df

def calculate_thresholds(thresholds_files):
    """Calculate classification thresholds from decoy datasets."""
    max_sps = float('-inf')
    min_sps = float('inf')
    for file in thresholds_files:
        df = pd.read_csv(file)
        df = calculate_sps(df)
        max_sps = max(max_sps, df['SPS'].max())
        min_sps = min(min_sps, df['SPS'].min())
    return {'L2': max_sps, 'L3': min_sps}

def run_cd_hit(fasta_input, fasta_output, threshold=0.8, word_size=5):
    """Run CD-HIT to cluster sequences."""
    cmd = f"cd-hit -i {fasta_input} -o {fasta_output} -c {threshold} -n { word_size}"
    subprocess.run(cmd, shell=True, check=True)

def align_with_hmm(fasta_file, hmm_file, output_sto_file):
    """Align sequences in a FASTA file to a profile HMM using hmmalign."""
    cmd = f"hmmalign --trim -o {output_sto_file} {hmm_file} {fasta_file}"
    run_external_command(cmd)

def convert_csv_to_fasta(csv_file, fasta_file):
    """Convert CSV file with sequence data to FASTA format, filtering for high-quality protein data."""
    df = pd.read_csv(csv_file)
    with open(fasta_file, 'w') as fasta:
        for _, row in df.iterrows():
            # Applying the three specified criteria
            if len(row['AA']) > 200 and row['AA'].startswith('M') and 'X' not in row['AA']: # TODO X in AA
                fasta.write(f">{row['name']}\n{row['AA']}\n")

def train_hmm(sto_file, hmm_output):
    """Train HMM model from Stockholm file."""
    cmd = f"hmmbuild {hmm_output} {sto_file}"
    run_external_command(cmd)

def combine_fasta_files(file_list, output_file):
    """Combine multiple FASTA files into a single FASTA file."""
    with open(output_file, 'w') as outfile:
        for fname in file_list:
            with open(fname) as infile:
                for line in infile:
                    outfile.write(line)

############################################################################################

def main_pipeline(thresholds_file, real_file, decoy_files, iterations=30, csv_dir = './input_csv', output_dir='./result_test', hmm_dir='./trained_hmm/', representative='./result_test/representative/'): 
    num_decoy = 10
    setup_directory(output_dir)
    setup_directory(csv_dir)
    setup_directory(hmm_dir)
    setup_directory(representative)

    previous_classifications = None # termination condition

    for iteration in range(1,iterations+1):
        thresholds_file = [f'./input_csv/i{iteration}/all_L7_decoy_iteration_{iteration}.fasta.csv']
        thresholds = calculate_thresholds(thresholds_file) # [1] cal SPS of i1_all_L7_decoys → get thresholds
        real_file = [f'./input_csv/i{iteration}/aligned_real.fasta.csv']

        for file in real_file: # for file in real_file + decoy_files:
            df = pd.read_csv(file)
            print(df)
            df = calculate_sps(df)
            df = classify_proteins(df, thresholds) # [2] Classify proteins
            processed_filename = os.path.join(output_dir, f'nonsplit_iteration_{iteration}_{os.path.basename(file)}')
            column_order = ['name', 'L2', 'L3', 'Length', 'SPS', 'Classification', 'AA']
            df = df[column_order]
            df.to_csv(processed_filename, index=False)
            
            # Split into L2-type and L3-type
            if file in real_file:  # Check if the current file is a real file to split
                l2_df = df[df['Classification'] == 'L2-type']
                l3_df = df[df['Classification'] == 'L3-type']
                non_df = df[df['Classification'] == 'Unclassified'] 

                # termination condition
                current_classifications = df[['name', 'Classification']].copy()
                if previous_classifications is not None:
                    # Check if the classifications have changed
                    if current_classifications.equals(previous_classifications):
                        print(f"No new classifications in iteration {iteration}. Terminating early.")
                        return  # Exit the function early if no new classifications are found

                # Update the previous_classifications with the current classifications
                previous_classifications = current_classifications

                ############################################################################################
                l2_filename = os.path.join(output_dir, f'L2_type_iteration_{iteration}_{os.path.basename(file)}')
                l3_filename = os.path.join(output_dir, f'L3_type_iteration_{iteration}_{os.path.basename(file)}')
                # all_L7_filename = os.path.join(output_dir, f'all_L7_iteration_{iteration}_{os.path.basename(file)}') 

                l2_df.to_csv(l2_filename, index=False)
                l3_df.to_csv(l3_filename, index=False)
                # all_L7_df = pd.concat([l2_df, l3_df], ignore_index=True)
                # all_L7_df.to_csv(all_L7_filename, index=False) # all_L7

                
                # Convert to FASTA and run CD-HIT
                l2_fasta = l2_filename.replace('.csv', '.fasta')
                l3_fasta = l3_filename.replace('.csv', '.fasta')
                # all_L7_fasta = all_L7_filename.replace('.csv', '.fasta')

                convert_csv_to_fasta(l2_filename, l2_fasta)
                convert_csv_to_fasta(l3_filename, l3_fasta)
                # convert_csv_to_fasta(all_L7_filename, all_L7_fasta)
              
                # Clustering with CD-HIT
                l2_clustered_fasta = l2_fasta.replace('.fasta', '_to_cluster.fasta')
                l3_clustered_fasta = l3_fasta.replace('.fasta', '_to_cluster.fasta')
                # all_L7_clustered = all_L7_fasta.replace('.fasta', '_to_cluster.fasta')
                
                run_cd_hit(l2_fasta, l2_clustered_fasta)
                run_cd_hit(l3_fasta, l3_clustered_fasta)
                # run_cd_hit(all_L7_fasta, all_L7_clustered) # all_L7

                # Align with previous iteration's HMM
                l2_hmm_file = os.path.join(hmm_dir, f'L2_i{iteration}.hmm')
                l3_hmm_file = os.path.join(hmm_dir, f'L3_i{iteration}.hmm')
                all_L7_hmm_file = os.path.join(hmm_dir, f'all_L7_i{iteration}.hmm') # all_L7
                


                # TODO: combine the two fasta l2_clustered l3_clustered be all_L7_clustered_fasta
                all_L7_clustered_fasta = os.path.join(output_dir, f'all_L7_iteration_{iteration}_to_cluster.fasta')
                combine_fasta_files([l2_clustered_fasta, l3_clustered_fasta], all_L7_clustered_fasta)

                l2_aligned_sto = l2_clustered_fasta.replace('_to_cluster.fasta', '_hmmalign_trim.sto')
                l3_aligned_sto = l3_clustered_fasta.replace('_to_cluster.fasta', '_hmmalign_trim.sto')
                all_L7_aligned_sto = all_L7_clustered_fasta.replace('_to_cluster.fasta', '_hmmalign_trim.sto')
            
                align_with_hmm(l2_clustered_fasta, l2_hmm_file, l2_aligned_sto)
                align_with_hmm(l3_clustered_fasta, l3_hmm_file, l3_aligned_sto)
                align_with_hmm(all_L7_clustered_fasta, all_L7_hmm_file, all_L7_aligned_sto) # all_L7
 
                ############################################################################################

                # Train HMM models
                l2_hmm_output = os.path.join(hmm_dir, f'L2_i{iteration+1}.hmm')
                l3_hmm_output = os.path.join(hmm_dir, f'L3_i{iteration+1}.hmm')
                all_L7_output = os.path.join(hmm_dir, f'all_L7_i{iteration+1}.hmm') 
                
                train_hmm(l2_aligned_sto, l2_hmm_output)
                train_hmm(l3_aligned_sto, l3_hmm_output)
                train_hmm(all_L7_aligned_sto, all_L7_output) 


                decoy_proteins_dir = os.path.join(output_dir, 'decoy_fasta')
                setup_directory(decoy_proteins_dir)
                destination = os.path.join(decoy_proteins_dir, f"all_L7_decoy_iteration_{iteration+1}.fasta")
                hmm_file = all_L7_output
                

                # 调用hmmemit命令
                command = f"hmmemit -o {destination} -N {num_decoy} {hmm_file}"
                subprocess.run(command, shell=True, check=True)

                ############################################################################################
                ############################################################################################
                ############################################################################################

        # 加載 L2 和 L3 模型
        with pyhmmer.plan7.HMMFile(l2_hmm_output) as hmm_file:
            hmm_L2 = next(hmm_file)

        with pyhmmer.plan7.HMMFile(l3_hmm_output) as hmm_file:
            hmm_L3 = next(hmm_file)

        # 創建字母表
        alphabet = pyhmmer.easel.Alphabet.amino()

        # decoy_fasta = glob.glob("./result_test/decoy_fasta/*.fasta")
        decoy_fasta = [f"./result_test/decoy_fasta/all_L7_decoy_iteration_{iteration+1}.fasta"]
        real_fasta = ["./real_fasta/aligned_real.fasta"]

        # csv_directory = os.path.join(output_dir, "csv")
        # setup_directory(csv_directory)  # 確保 csv 目錄已創建

        # 處理每個 FASTA 文件

        fasta_files = real_fasta + decoy_fasta
        for fasta_file in fasta_files:
            setup_directory(csv_dir+f"/i{iteration+1}")
            output_csv_path = os.path.join(csv_dir, f"i{iteration+1}/{os.path.basename(fasta_file)}.csv")

            with open(output_csv_path, 'w', newline='') as csvfile:
                fieldnames = ['name', 'L2', 'L3', 'Length', 'AA', 'Clade', 'Class']
                writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
                writer.writeheader()

                sequences = SeqIO.parse(fasta_file, "fasta")
                search_sequences = [TextSequence(name=bytes(seq.id, "utf-8"), sequence=str(seq.seq)) for seq in sequences]
                digital_search_sequences = [seq.digitize(alphabet) for seq in search_sequences]
                # sequences_block = DigitalSequenceBlock(alphabet, digital_search_sequences)

                # 將序列轉換為 DigitalSequence
                digital_search_sequences = [seq.digitize(alphabet) for seq in search_sequences]

                # 創建搜索管道並搜索 HMM
                background = pyhmmer.plan7.Background(alphabet)
                pipeline = pyhmmer.plan7.Pipeline(alphabet, background=background)

                sequences = pyhmmer.easel.DigitalSequenceBlock(alphabet, digital_search_sequences)

                # 使用 DigitalSequenceList 包裝序列
                hits_L2 = pipeline.search_hmm(query=hmm_L2, sequences=sequences)
                hits_L3 = pipeline.search_hmm(query=hmm_L3, sequences=sequences)

            # 構建名稱和分數字典以便於查找
                l2_scores = {hit.domains[0].alignment.target_name.decode(): hit.score for hit in hits_L2}
                l3_scores = {hit.domains[0].alignment.target_name.decode(): hit.score for hit in hits_L3}

                for seq in search_sequences:
                    name = seq.name.decode()
                    l2_score = l2_scores.get(name, "N/A")
                    l3_score = l3_scores.get(name, "N/A")
                    length, aa = len(seq.sequence), seq.sequence
                    
                    if fasta_file in real_fasta:
                        # Attempt to get the classification from the dataframe
                        classification2 = df.loc[df['name'] == name.split('|')[0], 'Classification'].values[0] if not df[df['name'] == name.split('|')[0]].empty else 'No class info available'
                        
                        writer.writerow({'name': name.split('|')[0], 'L2': f"{l2_score:.2f}" if l2_score != "N/A" else "N/A",
                                        'L3': f"{l3_score:.2f}" if l3_score != "N/A" else "N/A", 'Length': length, 'AA': aa,
                                        'Clade': name_to_lineage[name.split('|')[0]].split(',')[3], 'Class': classification2})
                    else:
                        writer.writerow({'name': name.split('-')[1], 'L2': f"{l2_score:.2f}" if l2_score != "N/A" else "N/A",
                                        'L3': f"{l3_score:.2f}" if l3_score != "N/A" else "N/A", 'Length': length, 'AA': aa,
                                        'Clade': 'Decoy', 'Class': 'Decoy'})


        print(f"########################## plot interaton {iteration} ##########################" ) 

        # 讀取並處理每個文件
        # 读取 CSV 文件
        file_path = f'./input_csv/i{iteration}/aligned_real.fasta.csv'  # 需要替换 {iteration} 为实际迭代数
        df = pd.read_csv(file_path)

        # 计算 x 和 y 值
        df['x'] = df['L2'] / df['Length']
        df['y'] = df['L3'] / df['Length']

        # 创建图表
        plt.figure(figsize=(10, 8))

        # 获取默认颜色循环
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colors = prop_cycle.by_key()['color']

        # 根据 Clade 分组并绘制
        for (clade, group), color in zip(df.groupby('Clade'), colors):
            plt.scatter(group['x'], group['y'], color=color, label=clade, alpha=0.6)

        plt.plot([0, 2], [0, 2], 'k--', label='Diagonal line')

        ############################################################################################
        simplified_classification = {k.split('|')[1]: v for k, v in classification.items()}
        simplified_L2_name = {name.split('|')[1] for name in L2_name}
        simplified_L3_name = {name.split('|')[1] for name in L3_name}
        # 颜色和标记设置
        colors_GT = {'L2': 'blue', 'L3': 'red'}
        markers_GT = {'Animal': 'o', 'Plant': 'x', 'Fungi': '^'}

        known_df = df[df['name'].apply(lambda x: simplified_classification.get(x, 'Unknown') != 'Unknown')]
    

        added_labels = set()

        for index, row in known_df.iterrows():
            if row['name'] in simplified_L2_name:
                group = 'L2'
                color_GT = colors_GT['L2']
            elif row['name'] in simplified_L3_name:
                group = 'L3'
                color_GT = colors_GT['L3']
            else:
                group = 'Unknown'  # 设置未知组别

            if group != 'Unknown':  # 只有当组别不是 Unknown 时才绘制
                category = simplified_classification.get(row['name'], 'Unknown')
                marker = markers_GT.get(category, '.')
                plot_label = f'{group} {category}'

                # 检查标签是否已经添加到图例中
                if plot_label not in added_labels:
                    plt.scatter(row['x'], row['y'], color=color_GT, marker=marker, label=plot_label)
                    added_labels.add(plot_label)
                else:
                    plt.scatter(row['x'], row['y'], color=color_GT, marker=marker)

        # 添加图例和标签
        plt.xlabel('L2 Score / Length')
        plt.ylabel('L3 Score / Length')
        plt.title(f'Iteration {iteration} with {num_decoy} decoys')
        plt.legend(title='Clade')
        plt.grid(True)
        plt.savefig(f"./png/i{iteration}.png")
        plt.close()

        # Additional plot: Proportion of Class within Clade
        if (iteration > 1):
            file_path = f'./input_csv/i{iteration}/aligned_real.fasta.csv'  # 需要替换 {iteration} 为实际迭代数
            df = pd.read_csv(file_path)

            # 计算 x 和 y 值
            df['x'] = df['L2'] / df['Length']
            df['y'] = df['L3'] / df['Length']

            # Additional plot: Proportion of Class within Clade
            plt.figure()
            filtered_df = df[df['Class'] != 'No class info available']
            class_clade_counts = filtered_df.groupby(['Clade', 'Class']).size().unstack(fill_value=0)
            class_clade_props = class_clade_counts.div(class_clade_counts.sum(axis=1), axis=0)

            class_clade_props.plot(kind='barh', stacked=True, color=['blue', 'red', 'green'])  # Use 'barh' for horizontal bars
            plt.title(f'Proportion of Clades within Classes (Iteration {iteration})')
            plt.ylabel('Clade')  # y-axis label for Class
            plt.xlabel('Proportion')  # x-axis label for Distribution
            plt.legend(title='Clade')
            plt.grid(True)
            plt.tight_layout()
            plt.savefig(f"./png/i{iteration}_class_proportion.png")
            plt.close()
    print("Pipeline completed successfully.")



# Run the main pipeline
thresholds_file = ['./input_csv/i1/all_L7_decoy.csv']
real_file = ['./input_csv/i1/aligned_real.csv']
decoy_files = ['./input_csv/i1/all_L7_decoy.csv', './input_csv/i1/L2_decoy.csv', './input_csv/i1/L3_decoy.csv']

import pandas as pd

# File paths
fasta_file_path = 'aligned_real.fasta'
tsv_file_path = 'LUC7p_species_phylogeny.tsv'

# Read the TSV file containing taxID and Named Lineage
tsv_df = pd.read_csv(tsv_file_path, sep='\t')
# Create a dictionary to map taxID to Named Lineage
taxid_to_lineage = {str(row['Taxid']): row['Named Lineage'] for index, row in tsv_df.iterrows()}

# Initialize a dictionary to hold the mapping from FASTA name to Named Lineage
name_to_lineage = {}

# Process the FASTA file
with open(fasta_file_path, 'r') as fasta_file:
    for line in fasta_file:
        line = line.strip()
        if line.startswith('>'):
            # Extract the identifier and remove the leading '>'
            identifier = line[1:]  # Remove '>' from the start of the header
            parts = identifier.split('|')
            name = parts[0]  # Assume the first part of the header is the name, now cleaned of '>'
            taxID = parts[-1].split(':')[-1]  # Extract taxID from the last part
            # Retrieve the Named Lineage using taxID
            named_lineage = taxid_to_lineage.get(taxID, 'NA,NA,NA,NA,NA')
            # Map the name to its corresponding Named Lineage
            name_to_lineage[name] = named_lineage




main_pipeline(thresholds_file, real_file, decoy_files)
