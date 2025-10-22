prior="L1"

# Create or clear the runtime log file
runtime_log="../../output/study2/signatureanalyzer_${prior}_runtime_2.csv"
echo "N,rep,seconds" > "$runtime_log"

for task_id in {1..20}
do
  bash study2_run_SignatureAnalyzer.sh $task_id $prior $runtime_log
done