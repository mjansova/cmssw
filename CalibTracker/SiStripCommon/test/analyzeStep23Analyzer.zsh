#!/bin/env zsh
#EPERADC X1 X2 X3
#declare -i m=1
#'s/MODE/'$(($i*20))'/g'
#

#for k in `seq 1 1` #7
#do
    #eadc=$((225+($k*5)))
    k=0
    eadc=247
    for i in `seq 1 20` #20
    do
        par2=$((0.04+($i*0.005)))
	for j in `seq 1 20` #20
	do
	    par3=$((0.001+(($j -1)*0.002)))
	    par1=$((1-(2*${par2})-(2*${par3})))
            sname="Coll2017_"${eadc}"_"${par2}"_"${par3}
            jobname="jobA_"${k}"_"${i}"_"${j}".sh"
	    #cat jobStep23Analyzer.sh | sed -e 's/FNAME/'$sname'/g' -e 's/EPERADC/'${eadc}'/g'  -e 's/X1/'${par1}'/g' > jobA_${k}_${i}_${j}.sh
	    cat jobStep23Analyzer.sh | sed -e 's/FNAME/'$sname'/g' -e 's/EPERADC/'${eadc}'/g' -e 's/X1/'${par1}'/g' -e 's/X2/'${par2}'/g' -e 's/X3/'${par3}'/g' -e 's/JOBN/'$jobname'/g'  > $jobname
            qsub -q cms_local_mdm jobA_${k}_${i}_${j}.sh 
	done
    done
#done
