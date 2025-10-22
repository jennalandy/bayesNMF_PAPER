# identify task id from command line arguments
task_id=$1
echo "task_id: $task_id"

# Read the CSV file and get the matrix value for the given task_id
matrix=$(awk -F',' -v id="$task_id" 'NR > 1 && $1 == id {gsub(/"/, "", $2); print $2}' pcawg_assignments.csv)
echo "matrix: $matrix"
name="${matrix%.csv}"
echo "name: $name"

# Construct the filename
filename="../../processed/PCAWG/matrices_nonhyper/${name}.csv"
output_dir="../../output/PCAWG/SignatureAnalyzer/${name}/"
mkdir -p "$output_dir"
echo "filename: $filename"

# Run the signatureanalyzer command with the constructed filename
signatureanalyzer -t matrix \
                  -n 10\
                  -o "$output_dir" \
                  --reference cosmic3 \
                  --objective poisson \
                  --prior_on_H L2 \
                  --prior_on_W L2 \
                  "$filename"