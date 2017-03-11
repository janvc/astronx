#!/bin/bash

echo "******************************************"
echo " running  003_wb"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_003_wb.inp
sleep 1

echo "******************************************"
echo " running  010_wb"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_010_wb.inp
sleep 1

echo "******************************************"
echo " running  020_wb"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_020_wb.inp
sleep 1

echo "******************************************"
echo " running  050_wb"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_050_wb.inp
sleep 1

echo "******************************************"
echo " running  050_ib"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_050_ib.inp
sleep 1

echo "******************************************"
echo " running  100_wb"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_100_wb.inp
sleep 1

echo "******************************************"
echo " running  100_ib"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_100_ib.inp
sleep 1

echo "******************************************"
echo " running  200_ib"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_200_ib.inp
sleep 1

echo "******************************************"
echo " running  500_ib"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_500_ib.inp
sleep 1

echo "******************************************"
echo " running  1k0_ib"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_1k0_ib.inp
sleep 1

echo "******************************************"
echo " running  2k0_ib"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_2k0_ib.inp
sleep 1

echo "******************************************"
echo " running  5k0_ib"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_5k0_ib.inp
sleep 1

echo "******************************************"
echo " running  10k_ib"
echo "******************************************"
time ../astronx_build/bin/astronx ../astronx/benchmark/bench_10k_ib.inp

exit 0

