# identify task id from command line arguments
task_id=$1
echo "task_id: $task_id"

prior=$2
echo "prior: $prior"

runtime_log=$3
echo "runtime_log: $runtime_log"

# Read the CSV file and get the matrix value for the given task_id
N=$(awk -F',' -v id="$task_id" 'NR > 1 && $1 == id {gsub(/"/, "", $2); print $2}' study2_sparse_assignments.csv)
echo "N: $N"

for rep in {1..10}; do
  name="N${N}_G64_rep${rep}"
  echo "name: $name"

  # Construct the filename
  filename="../../data/study2_sparse/${name}.csv"
  output_dir="../../output/study2_sparse/signatureanalyzer_${prior}/${name}/"

  # Check if the output directory already exists
  if [ -d "$output_dir" ]; then
    echo "Output directory $output_dir already exists. Skipping iteration."
    continue
  fi
  
  mkdir -p "$output_dir"
  echo "filename: $filename"

  # Run the signatureanalyzer command with the constructed filename
  start_time=$(date +%s)
  signatureanalyzer -t matrix \
                    -n 10\
                    -o "$output_dir" \
                    --reference cosmic3 \
                    --objective poisson \
                    --prior_on_H $prior \
                    --prior_on_W $prior \
                    "$filename"
  end_time=$(date +%s)
  runtime=$((end_time - start_time))

  # Append the runtime to the log file
  echo "$N,$rep,$runtime" >> "$runtime_log"
done