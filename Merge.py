import pandas as pd
import glob
import os

# ==============================================================================
# 1. ROBUST PATH SETUP
# ==============================================================================
# Get the root directory where this script is currently saved (e.g., the GBM folder)
script_dir = os.path.dirname(os.path.abspath(__file__))

print(f"📂 Script location: {script_dir}")
print(f"🔍 Searching for ML datasets across all subfolders...")

# ==============================================================================
# 2. FIND FILES (Recursive Search)
# ==============================================================================
# Search for any CSV file ending in "_ML_Ready.csv" in all subdirectories
search_pattern = os.path.join(script_dir, "**", "*_ML_Ready.csv")
all_files = glob.glob(search_pattern, recursive=True)

if not all_files:
    print("\n❌ NO ML_READY FILES FOUND!")
    print("---------------------------------------------------")
    print("1. Ensure you have run the R scripts for the datasets first.")
    print("2. Verify that the outputs (e.g., GSE61710_ML_Ready.csv) exist in their respective folders.")
    print("---------------------------------------------------")
    exit()

print(f"✅ Found {len(all_files)} dataset(s):")
for f in all_files:
    # Print just the folder and filename for cleaner output
    print(f"   - {os.path.relpath(f, script_dir)}")

# ==============================================================================
# 3. MERGE FILES
# ==============================================================================
dataframes = []

for filename in all_files:
    df = pd.read_csv(filename)
    # Ensure Sample_ID is a string to prevent merging issues
    df['Sample_ID'] = df['Sample_ID'].astype(str)
    dataframes.append(df)

print("\n🔄 Merging datasets (Finding common miRNAs)...")

# Join = 'inner' ensures we only keep miRNAs present in ALL files
master_df = pd.concat(dataframes, axis=0, join='inner', ignore_index=True)

# ==============================================================================
# 4. SAVE MASTER FILE
# ==============================================================================
# Save the output file in the root directory
output_file_path = os.path.join(script_dir, "MASTER_GBM_DATASET.csv")

# Reorder columns to ensure Sample_ID is first and Target is last
cols = list(master_df.columns)
if 'Sample_ID' in cols: cols.remove('Sample_ID')
if 'Target' in cols: cols.remove('Target')
final_order = ['Sample_ID'] + cols + ['Target']
master_df = master_df[final_order]

# Shuffle rows to prevent model bias based on ordered dataset ingestion
master_df = master_df.sample(frac=1, random_state=42).reset_index(drop=True)

master_df.to_csv(output_file_path, index=False)

print("\n=======================================================")
print(f"🎉 DONE! File saved to:\n{output_file_path}")
print(f"📊 Total Samples: {master_df.shape[0]}")
print(f"🧬 Common miRNAs: {master_df.shape[1] - 2}")
print("=======================================================")

if master_df.shape[1] - 2 < 5:
    print("⚠️  WARNING: Very few common miRNAs found.")
    print("    Check your naming conventions across datasets (e.g. hsa-miR-21 vs hsa-mir-21).")