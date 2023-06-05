#!/bin/bash

# Usage: ./window_survey job minW maxW
#
# First ./init a job with a dummy W.
# Then mv that folder to a new name without the W.
# If necessary, load/edit/save any other settings (eg. random seed).
# Then run this script.
#
MODULE="UniformWindows"

JOB=$1; shift
MINW=$1; shift
MAXW=$1; shift

SCRIPT="
    using $MODULE;
    minW=parse(Int,\"$MINW\");
    maxW=parse(Int,\"$MAXW\");
    for W in range(minW,maxW);
        load(\"$JOB\"; newdir=\"$JOB/W_\$W\", trace=false, state=false);
        _!.setup.init_Ω = _!.setup.ΩMAX;
        _!.setup.W = W;
        save();
    end;
"

julia --project=. -e "$SCRIPT"

for W in $(seq $MINW $MAXW)
do
    cp "$JOB/script.jl" "$JOB/W_$W"
    ./start "$JOB/W_$W"
done
