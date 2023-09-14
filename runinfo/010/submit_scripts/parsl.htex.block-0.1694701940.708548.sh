
export JOBNAME=$parsl.htex.block-0.1694701940.708548
set -e
export CORES=$(getconf _NPROCESSORS_ONLN)
[[ "1" == "1" ]] && echo "Found cores : $CORES"
WORKERCOUNT=1
FAILONANY=0
PIDS=""

CMD() {
process_worker_pool.py  --max_workers=9 -a phpmb48.local,10.154.7.254,152.7.255.206 -p 0 -c 1.0 -m None --poll 10 --task_port=54195 --result_port=54761 --logdir=/Users/jrabasc/Desktop/github/q2-quality-control/runinfo/010/htex --block_id=0 --hb_period=30  --hb_threshold=120 --cpu-affinity none --available-accelerators  --start-method spawn
}
for COUNT in $(seq 1 1 $WORKERCOUNT); do
    [[ "1" == "1" ]] && echo "Launching worker: $COUNT"
    CMD $COUNT &
    PIDS="$PIDS $!"
done

ALLFAILED=1
ANYFAILED=0
for PID in $PIDS ; do
    wait $PID
    if [ "$?" != "0" ]; then
        ANYFAILED=1
    else
        ALLFAILED=0
    fi
done

[[ "1" == "1" ]] && echo "All workers done"
if [ "$FAILONANY" == "1" ]; then
    exit $ANYFAILED
else
    exit $ALLFAILED
fi
