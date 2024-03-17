#!/bin/bash


exe="../4.1_friction/execute"
parfile="Parameters_LE.par"
basedir="../dataApril2015/gamma1e-6/"
mkdir -p $basedir
#basedir="/home/yffily/noAlign/data/5/2_noAlign_polydisperse_phasediag"
#inibasedir="/home/yffily/noAlign/data/5/1_noAlign_phasediag"
log="runMay2015N1024"

N_list="1024"
mu_list="0.1"
xi_list="0.1"
# phi_list="0.8275 0.83 0.8325 0.835"  
phi_list="0.8225 0.8250 0.8375 0.84" 
conf_list="1"
dt="2e-2"
polydis="0.2"
tsim="1e6"
tprint="1e+2"
strain_rate="1e-6"	

#######################################################################

## loop over parameters; create suitable directories; build dag file ##
for N in $N_list; do
	N_dir="$basedir/N_${N}"
	mkdir "$N_dir"
	for mu in $mu_list; do
		mu_dir="$N_dir/mu_${mu}"
		mkdir "$mu_dir"
		for xi in $xi_list; do
			xi_dir="$mu_dir/xi_${xi}"
			mkdir "$xi_dir"
			for phi in $phi_list; do
				phi_dir="$xi_dir/phi_${phi}"
				mkdir "$phi_dir"
				for k in $conf_list; do
					dir="$phi_dir/conf${k}"
					mkdir "$dir"
# 					 secho "tau=$tau, v=$v, phi=$phi, conf $k"
					echo "$exe" > "$dir/executable.log"
							
					newparfile="$dir/$parfile"
					cp "$parfile" "$newparfile"
				
					sed -i "s|outdir *= *[^\t# ]*|outdir = $dir|" "$newparfile"
					sed -i "s/np *= *[^\t# ]*/np = $N/" "$newparfile"
					sed -i "s/mu *= *[^\t# ]*/mu = $mu/" "$newparfile"
					sed -i "s/xivisc *= *[^\t# ]*/xivisc = $xi/" "$newparfile"
					sed -i "s/phi *= *[^\t# ]*/phi = $phi/" "$newparfile"
					sed -i "s|^dt *= *[^\t# ]*|dt = $dt|" "$newparfile"
					sed -i "s|^tsim *= *[^\t# ]*|tsim = $tsim|" "$newparfile"
					sed -i "s|^tprint *= *[^\t# ]*|tprint = $tprint|" "$newparfile"
					sed -i "s|^polydis *= *[^\t# ]*|polydis = $polydis|" "$newparfile"
					sed -i "s|^strain_rate *= *[^\t# ]*|strain_rate = $strain_rate|" "$newparfile"
					
					initype="Random"
					sed -i "s|initype *= *[^\t# ]*|initype = $initype|" "$newparfile"
					sed -i "s|iniangletype *= *[^\t# ]*|iniangletype = $initype|" "$newparfile"
# 					sed -i "s|seed_iniPos *= *[^\t# ]*|seed_iniPos = $[670*k+346]|" "$newparfile"
# 					sed -i "s|seed_iniRad *= *[^\t# ]*|seed_iniRad = $[354*k+384]|" "$newparfile"
					sed -i "s|seed_iniPos *= *[^\t# ]*|seed_iniPos = $[357*k+174]|" "$newparfile"
					sed -i "s|seed_iniRad *= *[^\t# ]*|seed_iniRad = $[489*k+159]|" "$newparfile"

					echo ${exe} ${newparfile} 
					${exe} ${newparfile} &
				done
			done
		done
	done
done