# TDE _ VORONOI_analysis_LAMMPS
Evaluation of the threshold displacment energy for Tungsten. 
# displacement simulation, lammps input# 
variable proj1  string       PKA_D0000_E080_T0000
variable pkaid1 equal        780
variable pkavx1 equal        15.752654
variable pkavy1 equal        -27.281989
variable pkavz1 equal        -81.841847
variable res_file1 string    /home/mosab/resources/restart.W-30K.100000
variable pot_file1 string    /home/mosab/resources/W_BN.eam.fs
variable atom_symbol1 string W
variable atom_mass1   equal  183.840000
variable tmax1 equal         0.002000
variable xmax1 equal         0.010000
variable nrun1 equal         7000
variable res_freq1 equal     1000
variable thermo_freq1 equal  5
variable dt_freq1     equal  10
variable discard_nrun1  equal 0
variable discard_tstep1 equal 0.001000


units metal
boundary p p p
atom_style atomic
atom_modify map array
neighbor 1.0 bin                # sets maximum neighbor search radius to cutoff+value, using bin-sort algorithm
neigh_modify delay 5 check yes  # checks if neighbor list should be rebuilt every 5 steps
read_restart ${res_file1}
compute 5 all voronoi/atom occupation
pair_style eam/fs
pair_coeff * * ${pot_file1} ${atom_symbol1}
mass 1 ${atom_mass1}
group PKA id  == ${pkaid1} 
fix 1 all nve
timestep ${discard_tstep1}    # to change the timing of the initioation of the recoil event
run      ${discard_nrun1}     # to change the timing of the initioation of the recoil event
reset_timestep 0
thermo  ${thermo_freq1}
log     log.${proj1}
thermo_style custom step time  dt time temp ke  pe etotal press
timestep  ${tmax1}
run 100

reset_timestep 0
dump 1  PKA custom 1000 dump.PKA_${proj1} id type x y z xs ys zs vx vy vz
dump 2  all custom 1000  dump.voronoi_${proj1} c_5[1]
write_data data.bofore_event
#restart ${res_freq1} restart.${proj1}
velocity PKA set  ${pkavx1} ${pkavy1} ${pkavz1}   units box 
fix 2 all dt/reset ${dt_freq1} NULL  ${tmax1} ${xmax1} units box #fix ID group-ID dt/reset N Tmin Tmax Xmax keyword values ...
fix 1com all recenter INIT INIT INIT
#dump 3  all  xyz 5000   dump.trajectories_${proj1} #type x y z xs ys zs vx vy vz
dump 3  PKA custom 1000 dump.PKA_d_${proj1} id type x y z xs ys zs vx vy vz
dump 4 all custom  100 trajectory.dump.${proj1} id type  x y z ix iy iz vx vy vz     
run ${nrun1}
write_data data.after_event_20ps.${pka_ene}
