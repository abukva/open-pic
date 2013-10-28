for f in parametri/*.conf
do
	mpirun -np $1 ./simulacija $f
done
