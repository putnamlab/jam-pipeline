#!/usr/bin/env bash

for f in kmerpipe.py fastq2pfa.pl gzcat gzip # testmissing
do
    if [ `which $f` ] ; then echo "$f 	-- OK."; else echo "$f 	-- NOT FOUND."; fi
done

echo "If things weren't found, check if " ${0/check_env.bash/} " is in your PATH"

