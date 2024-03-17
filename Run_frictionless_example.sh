#!/bin/bash


exe="../4.1_friction/execute"
parfile="Parameters_LE_Frictionless.par"
basedir="../DataFricless/"
mkdir -p $basedir

N_list="64 256 1024"
xi_list="1"
gam_list="1e-6"
phi_list="0.835 0.8375 0.84 0.841 0.842 0.843 0.845 0.85" # 
Nconf=1
dt="2e-2"
polydis="0.2"
tsim="1e6"

u=5873762
## loop over parameters; create suitable directories; build dag file ##
for N in $N_list; do
	N_dir="$basedir/N_${N}"
	mkdir "$N_dir"
	for gam in $gam_list; do
		gam_dir="$N_dir/gam_${gam}"
		mkdir "$gam_dir"
		for xi in $xi_list; do
			xi_dir="$gam_dir/xi_${xi}"
			mkdir "$xi_dir"
			for phi in $phi_list; do
				phi_dir="$xi_dir/phi_${phi}"
				mkdir "$phi_dir"
				for k in $(seq $Nconf); do
					dir="$phi_dir/conf${k}"
					mkdir "$dir"
					echo "$dir"
					echo "tau=$tau, v=$v, phi=$phi, conf $k"
							
					newparfile="$dir/$parfile"
					logfile="$dir/log.dat"
					cp "$parfile" "$newparfile"
				
					sed -i "s|outdir *= *[^\t# ]*|outdir = $dir|" "$newparfile"
					sed -i "s/np *= *[^\t# ]*/np = $N/" "$newparfile"
					sed -i "s/strain_rate *= *[^\t# ]*/strain_rate = $gam/" "$newparfile"
					sed -i "s/xivisc *= *[^\t# ]*/xivisc = $xi/" "$newparfile"
					sed -i "s/phi *= *[^\t# ]*/phi = $phi/" "$newparfile"
					sed -i "s|^dt *= *[^\t# ]*|dt = $dt|" "$newparfile"
					sed -i "s|^tsim *= *[^\t# ]*|tsim = $tsim|" "$newparfile"
					sed -i "s|^polydis *= *[^\t# ]*|polydis = $polydis|" "$newparfile"			
					
					initype="Random"
					sed -i "s|initype *= *[^\t# ]*|initype = $initype|" "$newparfile"
					sed -i "s|iniangletype *= *[^\t# ]*|iniangletype = $initype|" "$newparfile"
					sed -i "s|seed_iniPos *= *[^\t# ]*|seed_iniPos = $u|" "$newparfile"
					sed -i "s|seed_iniRad *= *[^\t# ]*|seed_iniRad = $[u+1]|" "$newparfile"
					echo $exe $newparfile > $logfile
					$exe $newparfile > $logfile &
					
					u=$[$u+2]
				done
			done
		done
	done
done

