#!/bin/bash
out='out_KS.txt' 
 rm kin_*
 rm *_kin.dat
 rm *_pot.dat
 rm A*.dat
 rm C_*.dat
 rm H*.dat
 rm O*.dat
 rm N*.dat
 grep "| Total number of basis functions :" $out |tee for_stm.txt
 grep "| Number of Kohn-Sham states (occupied + empty):" $out
 grep -n "State" $out | tail -1
 echo "Now wrile commmand: "
 echo "sed -n start_line,end_line p output.file > all_states.dat"
 echo "where end_line = startline + number of states + 2"


