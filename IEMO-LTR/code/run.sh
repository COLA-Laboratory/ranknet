#!/bin/bash
beginTime=`date +%s`

mkdir ../log

start_index=0
run=3
elicitationInterval=10
inquiriesNum=10
stepSize=0.2
noise=-1
thread=$run
save=10

config=configFiles/test_problem/dtlz/dtlz2_3d.txt
beforeTime=`date +%s`
./EMOC -file $config -algorithm moead_ltr -W 1,1,1 -save $save -start_index $start_index -run $run -thread $thread -elicitationInterval $elicitationInterval -inquiriesNum $inquiriesNum -stepSize $stepSize -noise $noise >../log/log.$beforeTime \
&&
echo "./EMOC -file $config -algorithm moead_ltr -W 1,1,1 -save $save -start_index $start_index -run $run -thread $thread -elicitationInterval $elicitationInterval -inquiriesNum $inquiriesNum -stepSize $stepSize -noise $noise >../log/log.$beforeTime"
afterTime=`date +%s`
echo "finished in $(($afterTime-$beforeTime)) seconds."

wait
endTime=`date +%s`
echo "all tasks run in $(($endTime-$beginTime)) seconds."
exit