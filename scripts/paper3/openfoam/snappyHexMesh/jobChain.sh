#!/bin/bash
jid1=$(sbatch --parsable step1.sh)
jid2=$(sbatch --parsable --dependency=afterok:$jid1 step3.sh)
jid3=$(sbatch --parsable --dependency=afterok:$jid2 reRun.sh)
jid4=$(sbatch --parsable --dependency=afterok:$jid3 reRun.sh)
echo "Submitted: $jid1 -> $jid2 -> $jid3 -> $jid4"
