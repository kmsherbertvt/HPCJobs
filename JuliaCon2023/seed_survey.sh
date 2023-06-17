#!/bin/bash

# Usage: ./inspect jobdir1 jobdir2 ...
MODULE="JuliaCon2023"

JOB=$1; shift
MINSEED=$1; shift
MAXSEED=$1; shift

SCRIPT="
    using $MODULE;
    minseed=parse(Int,\"$MINSEED\");
    maxseed=parse(Int,\"$MAXSEED\");
    for seed in range(minseed,maxseed);
        load(\"$JOB\"; newdir=\"$JOB/seed_\$seed\", trace=false, state=false);
        _!.setup.seed = seed;
        save();
    end;
"

julia --project=. -e "$SCRIPT" $SEED

for SEED in $(seq $MINSEED $MAXSEED)
do
    cp "$JOB/script.jl" "$JOB/seed_$SEED"
    ./start "$JOB/seed_$SEED"
done
