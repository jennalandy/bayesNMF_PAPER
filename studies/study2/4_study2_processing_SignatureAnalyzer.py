import pandas as pd
import os

# process L1
results_dir = "../../output/study2/signatureanalyzer_L1"
results_list = os.listdir(results_dir)
for results in results_list:
  file_path = os.path.join(results_dir, results, 'nmf_output.h5')
  try:
    print(file_path)
    W = pd.read_hdf(file_path, 'W')
    W.to_csv(os.path.join(results_dir, results, "W.csv"))
    H = pd.read_hdf(file_path, 'H')
    H.to_csv(os.path.join(results_dir, results, "H.csv"))
  except Exception as e:
    print(file_path, 'can\'t open')

# process L2
results_dir = "../../output/study2/signatureanalyzer_L2"
results_list = os.listdir(results_dir)
for results in results_list:
  file_path = os.path.join(results_dir, results, 'nmf_output.h5')
  csv_path = os.path.join(results_dir, results, "W.csv")
  if os.path.exists(csv_path):
    continue
  try:
    print(file_path)
    W = pd.read_hdf(file_path, 'W')
    W.to_csv(os.path.join(results_dir, results, "W.csv"))
    H = pd.read_hdf(file_path, 'H')
    H.to_csv(os.path.join(results_dir, results, "H.csv"))
  except Exception as e:
    print(file_path, 'can\'t open')