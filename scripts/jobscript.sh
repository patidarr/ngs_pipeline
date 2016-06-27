#!/bin/bash
# properties = {properties}
set -x

# propagate the TERM signal to snakemake
trap 'kill -TERM $PID; touch {jobfailed}' SIGTERM

{exec_job} &
PID=$!
# wait exits immediately if the script is sent sigterm
# therefore, there has to be another wait in the trap
wait $PID
trap - SIGTERM
wait $PID
exit $?
