download_data_dir=/net/o3/hymet_nobackup/heimc/data/pgw/download
subdom_data_dir=/net/o3/hymet_nobackup/heimc/data/pgw/download/subdomain

box=-73,37,-42,34

# Note: use ps,hurs,tas from Amon and zg from Emon

var_names=(ta hur ua va)
experiments=(historical ssp585)

for var_name in ${var_names[@]}; do
    echo $var_name
    for exp in ${experiments[@]}; do
        echo $exp
        if [[ "$exp" == "historical" ]]; then
            years=(19850101-19891231 19900101-19941231 \
                   19950101-19991231 20000101-20041231 \
                   20050101-20091231 20100101-20141231)
            years=(19850101-19891231 19900101-19941231)
        elif [[ "$exp" == "ssp585" ]]; then
            years=(20700101-20751231 20750101-20791231 \
                   20800101-20841231 20850101-20891231 \
                   20900101-20941231 20950101-20991231)
            years=(20700101-20751231 20750101-20791231)
        fi
        for year in ${years[@]}; do
            echo $year
            cdo sellonlatbox,$box \
                $download_data_dir/${var_name}_CFday_MPI-ESM1-2-HR_${exp}_r1i1p1f1_gn_${year}.nc \
                $subdom_data_dir/${var_name}_CFday_MPI-ESM1-2-HR_${exp}_r1i1p1f1_gn_${year}.nc
        done
    done
done
