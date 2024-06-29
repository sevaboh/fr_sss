g++ -O3 fr_soil_strain_stress.cpp -o fr_sss
#g++ -fopenmp -O3 fr_soil_strain_stress.cpp -o fr_sss
order=1.0
order2=1.0
biot=1.0
biot2=1.0
#time ./fr_sss testing_mode 1 NB 50 debug_level 0 b2 $biot2 biot $biot beta $order2 alpha $order Zoutstep 1 Tm 864000 Om 1 LB 1 pick_eps 1e-12 ls_eps 1e-12 ls_e2 1e-14 fr_eps 1e-20 ls_max_iter 1200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt max_tau 400 st_tau 0.01 ls_min_tau 0.00000001
#exit
time ./fr_sss testing_mode 0 NB 150 debug_level 1 is_irr 1 b2 $biot2 biot $biot beta $order alpha $order Zoutstep 1 Tm 30000000 Om 1 LB 1 pick_move 0.5 pick_eps 1e-10 ls_eps 1e-12 ls_e2 1e-14 fr_eps 1e-20 ls_max_iter 4000 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt max_tau 3600 st_tau 0.01 ls_min_tau 0.00000001
exit
time ./fr_sss testing_mode 1 NB 25 debug_level 0 b2 $biot2 biot $biot beta $order2 alpha $order Zoutstep 1 Tm 864000 Om 3600 LB 1 pick_eps 1e-12 ls_eps 1e-12 ls_e2 1e-14 fr_eps 1e-20 ls_max_iter 1200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt max_tau 800 st_tau 0.01 ls_min_tau 0.00000001 > 25.txt
time ./fr_sss testing_mode 1 NB 50 debug_level 0 b2 $biot2 biot $biot beta $order2 alpha $order Zoutstep 1 Tm 864000 Om 3600 LB 1 pick_eps 1e-12 ls_eps 1e-12 ls_e2 1e-14 fr_eps 1e-20 ls_max_iter 1200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt max_tau 400 st_tau 0.01 ls_min_tau 0.00000001 > 50.txt
time ./fr_sss testing_mode 1 NB 100 debug_level 0 b2 $biot2 biot $biot beta $order2 alpha $order Zoutstep 1 Tm 864000 Om 3600 LB 1 pick_eps 1e-12 ls_eps 1e-12 ls_e2 1e-14 fr_eps 1e-20 ls_max_iter 1200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt max_tau 200 st_tau 0.01 ls_min_tau 0.00000001 > 100.txt
time ./fr_sss testing_mode 1 NB 150 debug_level 0 b2 $biot2 biot $biot beta $order2 alpha $order Zoutstep 1 Tm 864000 Om 3600 LB 1 pick_eps 1e-12 ls_eps 1e-12 ls_e2 1e-14 fr_eps 1e-20 ls_max_iter 1200 vgm vgm1.txt et_file ET2.txt prec_file precipitation2.txt gw_file gw_depth3.txt max_tau 133 st_tau 0.01 ls_min_tau 0.00000001 > 150.txt
