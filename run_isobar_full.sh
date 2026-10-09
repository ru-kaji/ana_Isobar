#!/bin/bash

trap "kill 0" EXIT

source /usr/local/root/root_v6.22.02/bin/thisroot.sh

# 1. 解析ループ（並列実行）
file_list=""
for i in {3..8}
do
    runnum="1900002${i}"
    output_file="OutAnaTree_${runnum}.root"
    
    echo "Starting Analysis for Run ${runnum}..."
    ./Run "${runnum}.lst" "${runnum}" 2018 > "log_${runnum}.txt" 2>&1 &
    
    file_list="$file_list $output_file"
done

echo "Waiting for all analysis processes to finish..."
wait

# 2. 合算 (hadd)
combined_out="OutAnaTree_isobar_mycut_R1p0_full_phiminus.root"
echo "Merging files into $combined_out..."
hadd -f "$combined_out" $file_list

# 3. 合算データに対して描画を実行
echo "Generating final plots from combined data..."

echo "Done! All processes finished."

trap - EXIT

echo "run_isobar_full.sh 23-26" | mail -s "Finished" s2520168@u.tsukuba.ac.jp

sleep 7
# read dummy
