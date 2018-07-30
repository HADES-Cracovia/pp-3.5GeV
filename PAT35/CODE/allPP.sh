filelist=$1
INPUT_DIR="/lustre/hades/user/beatam/dst_simFiz/files/pp_elastic";
OUPUT_DIR="/lustre/hades/user/przygoda/FILES/PP/SIM/pp_elastic/";
for item in $(cat $filelist)
do
    cp $INPUT_DIR/${item}.root /tmp
    /u/przygoda/PAT/anapp ${item}.root
    mv /tmp/${item}_hadron_out.root $OUPUT_DIR
    rm -rf /tmp/${item}*.root
done


