#!/bin/bash

cd selectedSeismograms
for ievt in `seq 59 1 113`;do
   echo $ievt
   cd $ievt
   mcpf=`ls 2*.mcp`
   mcpf_z=`echo $mcpf |awk -F. '{printf"%s.%s.BHZ.mcp",$1,$2}'`
   mv $mcpf $mcpf_z
   cd ..
done
