#!/bin/bash

#exc="addqueue -q cmb -s -n 1x12 -m 2.0 /usr/local/shared/python/2.7.6-gcc/bin/python"
exc="python"
nsims=1000

for nside in 64 256
do
    for msk in 1
    do
	for cont in 0 1
	do
	    for aposize in 1.
	    do
		command="${exc} check_sph.py ${nside} ${msk} ${cont} ${nsims} 0 ${aposize}"
		echo ${command}
		${command}
	    done
	done
    done
done

for msk in 1
do
    for cont in 0 1
    do
	for aposize in 0. 0.1
	do
	    command="${exc} check_flat.py ${msk} ${cont} ${nsims} 0 ${aposize}"
	    echo ${command}
	    #${command}
	done
    done
done

exc="python"
for nside in 64 256
do
    for msk in 1
    do
	for cont in 0 1
	do
	    for aposize in 1.
	    do
		command="${exc} check_sph.py ${nside} ${msk} ${cont} ${nsims} 1 ${aposize}"
		echo ${command}
		${command}
	    done
	done
    done
done

for msk in 1
do
    for cont in 0 1
    do
	for aposize in 0. 0.1
	do
	    command="${exc} check_flat.py ${msk} ${cont} ${nsims} 1 ${aposize}"
	    echo ${command}
	    #${command}
	done
    done
done
