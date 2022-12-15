CHECKER=$1
bins_dir=$2
REPORT=${bins_dir}_validate_report.txt
rm -rf ${REPORT}
touch ${REPORT}
# CHECKER=/home/mabuelanin/dib-dev/kSpider_bins/build/check_bin
no_bins=$(find ${bins_dir} -printf \\n | wc -l)
COUNTER=1
for bin in ${bins_dir}/*.bin;
do
    let COUNTER++;
    echo "${COUNTER}/${no_bins}";
    result=$(${CHECKER} ${bin} 2>&1);
    if [[ "${result}" == *"VALID_BIN"* ]]; 
    then
        echo -e "${bin} | ${result}" >> ${REPORT}
    else
        echo "${bin} | INVALID" >> ${REPORT}
    fi
done