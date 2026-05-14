g++ -O3 -I/usr/include/python3.10 -DNDEBUG -fstack-protector-strong -fwrapv fr_soil_strain_stress_particle_transport.cpp -o fr_sss_pt -lcrypt -ldl -lm -lpython3.10 
g++ -O3 -fopenmp -I/usr/include/python3.10 -DNDEBUG -fstack-protector-strong -fwrapv fr_soil_strain_stress_particle_transport.cpp -o fr_sss_pt_omp -lcrypt -ldl -lm -lpython3.10 
order=0.9
order2=0.9
biot=1.0
biot2=1.0
time ./fr_sss_pt_omp NB 150 debug_level 0 is_irr 1 b2 $biot2 biot $biot beta $order alpha $order Zoutstep 1 Tm 30000000 Om 3600 LB 1 pick_move 0.5 pick_eps 1e-10 ls_eps 1e-12 ls_e2 1e-14 fr_eps 1e-20 ls_max_iter 4000 vgm vgm2.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt max_tau 3600 st_tau 0.01 ls_min_tau 0.00000001 porosity 0.44 > out_pt_omp.txt
mv out_H.txt out_Hptomp.txt
mv out_V.txt out_Vptomp.txt
mv out_Vz.txt out_Vzptomp.txt
mv out_Z.txt out_Zptomp.txt
mv out_n.txt out_nptomp.txt
mv out_P_0.txt out_P0omp.txt
mv out_P_1.txt out_P1omp.txt
mv out_P_2.txt out_P2omp.txt
mv out_vgm.txt out_vgmomp.txt
time ./fr_sss_pt NB 150 debug_level 0 is_irr 1 b2 $biot2 biot $biot beta $order alpha $order Zoutstep 1 Tm 30000000 Om 3600 LB 1 pick_move 0.5 pick_eps 1e-10 ls_eps 1e-12 ls_e2 1e-14 fr_eps 1e-20 ls_max_iter 4000 vgm vgm2.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt max_tau 3600 st_tau 0.01 ls_min_tau 0.00000001 porosity 0.44 > out_pt.txt
mv out_H.txt out_Hpt.txt
mv out_V.txt out_Vpt.txt
mv out_Vz.txt out_Vzpt.txt
mv out_Z.txt out_Zpt.txt
mv out_n.txt out_npt.txt

g++ -O3 fr_soil_strain_stress_v2.cpp -o fr_sss
g++ -O3 -fopenmp fr_soil_strain_stress_v2.cpp -o fr_sss_omp
time ./fr_sss_omp testing_mode 0 NB 150 debug_level 0 is_irr 1 b2 $biot2 biot $biot beta $order alpha $order Zoutstep 1 Tm 30000000 Om 3600 LB 1 pick_move 0.5 pick_eps 1e-10 ls_eps 1e-12 ls_e2 1e-14 fr_eps 1e-20 ls_max_iter 4000 vgm vgm2.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt max_tau 3600 st_tau 0.01 ls_min_tau 0.00000001 porosity 0.44 rho_s 1287 > out_omp.txt 
mv out_H.txt out_Homp.txt
mv out_V.txt out_Vomp.txt
mv out_Vz.txt out_Vzomp.txt
mv out_Z.txt out_Zomp.txt
mv out_n.txt out_nomp.txt

time ./fr_sss testing_mode 0 NB 150 debug_level 0 is_irr 1 b2 $biot2 biot $biot beta $order alpha $order Zoutstep 1 Tm 30000000 Om 3600 LB 1 pick_move 0.5 pick_eps 1e-10 ls_eps 1e-12 ls_e2 1e-14 fr_eps 1e-20 ls_max_iter 4000 vgm vgm2.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt max_tau 3600 st_tau 0.01 ls_min_tau 0.00000001 porosity 0.44 rho_s 1287 > out.txt 
