#!/bin/env zsh
#EPERADC X1 X2 X3
#declare -i m=1
#'s/MODE/'$(($i*20))'/g'
#

#for k in `seq 1 7` #7
#do
    #eadc=$((225+($k*5)))
    eadc=247
    k=0
    for i in `seq 1 20` #20
    do
        par2=$((0.04+($i*0.005)))
	for j in `seq 1 20` #20
	do
	    par3=$((0.001+(($j -1)*0.002)))
	    par1=$((1-(2*${par2})-(2*${par3})))
            sname="Coll2017_"${eadc}"_"${par2}"_"${par3}
            smod="TOB"
            ldown=1
            lup=4
            side="mono"
            #side="stereo"
	    echo | cat jobChi2.sh | sed -e 's/FNAME/'$sname'/g' -e 's/EPERADC/'${eadc}'/g' -e 's/X1/'${par1}'/g' -e 's/X2/'${par2}'/g' -e 's/X3/'${par3}'/g' -e 's/MOD/'$smod'/g'  -e 's/LDOWN/'${ldown}'/g' -e 's/LUP/'${lup}'/g' -e 's/SIDE/'$side'/g'
	done
    done
#done


#if [ $j -eq 3 ]
#then
 #   m=4
#elif [ $j -eq 4 ]
#then
#    m=3 
#elif [ $j -eq 5 ]
#then
#    m=6 
#elif [ $j -eq 6 ]
#then
#    m=7
#else
#    echo Non defined partition!
#fi 
#
#    for k in `seq 1 $m`
#    do
#    for i in `seq 1 50` ; do cat jobAnalysis.sh | sed -e 's/MODB/'$(((20*($i-1))+1))'/g' -e 's/MODE/'$(($i*20))'/g'  > jobA.sh; echo qsub -q cms_local_mdm  jobA.sh ; done
#    done
#
#done
