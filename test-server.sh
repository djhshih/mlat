path=src

host=127.0.0.1
port=6667
ntries=5


$path/blatd || true
$path/blatc || true


$path/blatd start $host $port data/hg38_tp53.2bit &
blatd_pid=$!
echo "Started server at $host:$port with pid $blatd_pid"

for r in {1..$ntries}; do
	sleep 1
	if $path/blatc $host $port . data/reads_tp53.fa out/out-1.1.c.psl; then
		# query succeeded
		break
	fi
done

kill $blatd_pid
echo "Shutdown server at $host:$port with pid $blatd_pid"

