prior="L2"

# Create or clear the runtime log file
runtime_log="../../output/study2_sparse/signatureanalyzer_${prior}_runtime_4.csv"
echo "N,rep,seconds" > "$runtime_log"

for task_id in {1..20}
do
  bash study2_sparse_run_SignatureAnalyzer.sh $task_id $prior $runtime_log
done