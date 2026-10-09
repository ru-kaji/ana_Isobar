#!/bin/bash

trap "kill 0" EXIT

source /usr/local/root/root_v6.22.02/bin/thisroot.sh

# 1. 解析ループ（並列実行）
file_list=""
for i in {1..2}
do
    runnum="1900000${i}"
    output_file="OutAnaTree_${runnum}.root"
    
    echo "Starting Analysis for Run ${runnum}..."
    ./Run "${runnum}.lst" "${runnum}" 2018 > "log_${runnum}.txt" 2>&1 &
    
    file_list="$file_list $output_file"
done

echo "Waiting for all analysis processes to finish..."
wait

# 2. 合算 (hadd)
combined_out="OutAnaTree_isobar_2root.root"
echo "Merging files into $combined_out..."
hadd -f "$combined_out" $file_list

# 3. 合算データに対して描画を実行
echo "Generating final plots from combined data..."

echo "Done! All processes finished."

trap - EXIT

sleep 7
# read dummy
