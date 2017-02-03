#!/bin/bash

path=src

host=127.0.0.1
port=6667
ntries=5

mkdir -p out/server


$path/blatd || true
$path/blatc || true


$path/blatd start $host $port data/hg38_tp53.2bit &
blatd_pid=$!
echo "Started server at $host:$port with pid $blatd_pid"

for r in {1..$ntries}; do
	sleep 1
	if $path/blatc $host $port . data/reads_tp53.fa out/server/out-1.1.psl; then
		# query succeeded
		break
	fi
done

kill $blatd_pid
echo "Shutdown server at $host:$port with pid $blatd_pid"

echo
echo "Checking outputs..."
echo

cases=( \
	1.1.psl \
)

for x in ${cases[@]}; do
	echo "diff out/server/out-${x} test/server/test-${x}"
	diff out/server/out-${x} test/server/test-${x}
done

