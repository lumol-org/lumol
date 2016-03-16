#!/bin/bash

set -e

if [ ! -f Cargo.toml ]
then
    echo "Error: $0 must be run from the crate root directory"
fi

multirust override nightly 2> /dev/null


DATETIME=$(git log -n1 --format="%ad" --date="format:%F-%R")
RUSTC_INFO=$(rustc -vV)
TARGET=$(echo $RUSTC_INFO | perl -ne '/host: ([\w\d_-]*)/ && print $1')
OUTPUT="benches/results/$DATETIME-$TARGET.bench"

# Build the needed binaries
cargo bench --no-run

echo "#" $(rustc -V) > $OUTPUT
echo "#" $(git log -n1 --format="%h %ai %s") >> $OUTPUT
# Only run the benchmarcks in `benches`
for file in benches/*.rs
do
    name=$(basename $file | cut -f1 -d '.')
    cargo bench --bench $name | tee -a $OUTPUT
done
echo "Results written to $OUTPUT"

multirust remove-override
