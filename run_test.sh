#!/bin/bash

echo "tt"
./NaMaster -map test/map_t.fits -map_2 test/map_t.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 0 -pol_2 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_tt_scal -nlb 4
echo "te"
./NaMaster -map test/map_t.fits -map_2 test/map_e.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 0 -pol_2 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_te_scal -nlb 4
echo "tb"
./NaMaster -map test/map_t.fits -map_2 test/map_b.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 0 -pol_2 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_tb_scal -nlb 4
echo "et"
./NaMaster -map test/map_e.fits -map_2 test/map_t.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 0 -pol_2 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_et_scal -nlb 4
echo "ee"
./NaMaster -map test/map_e.fits -map_2 test/map_e.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 0 -pol_2 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_ee_scal -nlb 4
echo "eb"
./NaMaster -map test/map_e.fits -map_2 test/map_b.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 0 -pol_2 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_eb_scal -nlb 4
echo "bt"
./NaMaster -map test/map_b.fits -map_2 test/map_t.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 0 -pol_2 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_bt_scal -nlb 4
echo "be"
./NaMaster -map test/map_b.fits -map_2 test/map_e.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 0 -pol_2 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_be_scal -nlb 4
echo "bb"
./NaMaster -map test/map_b.fits -map_2 test/map_b.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 0 -pol_2 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_bb_scal -nlb 4
echo "tp"
./NaMaster -map test/map_t.fits -map_2 test/map_qu.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 0 -pol_2 1 -cl_noise test/cl_noise_2.dat -coupling none -out test/cl_tp_spin -nlb 4
echo "pp"
./NaMaster -map test/map_qu.fits -map_2 test/map_qu.fits -mask test/mask.fits -mask_2 test/mask.fits -pol 1 -pol_2 1 -cl_noise test/cl_noise_4.dat -coupling none -out test/cl_pp_spin -nlb 4
echo "tt1"
./NaMaster -map test/map_t.fits -mask test/mask.fits -pol 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_tt_scal_1 -nlb 4
echo "ee1"
./NaMaster -map test/map_e.fits -mask test/mask.fits -pol 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_ee_scal_1 -nlb 4
echo "bb1"
./NaMaster -map test/map_b.fits -mask test/mask.fits -pol 0 -cl_noise test/cl_noise_1.dat -coupling none -out test/cl_bb_scal_1 -nlb 4
echo "pp1"
./NaMaster -map test/map_qu.fits -mask test/mask.fits -pol 1 -cl_noise test/cl_noise_4.dat -coupling none -out test/cl_pp_spin_1 -nlb 4
