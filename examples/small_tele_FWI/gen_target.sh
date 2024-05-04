cat > grid << eof
833970.438 -44272.6367 -80000
2230.7389111111406 2216.021092499999 2000
91 41 41
eof

for name in vs vp rho; do
    mpirun -np 4 ../../bin/xmodel_gll grid ${name} ./${name}_target.npz model_target/
done

rm grid