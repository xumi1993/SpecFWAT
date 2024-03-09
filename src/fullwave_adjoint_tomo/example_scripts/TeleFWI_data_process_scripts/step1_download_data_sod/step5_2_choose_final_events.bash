#!/bin/bash
#cl=d0-30  # catalog <30km or >200km
#cl=d200-10000  # catalog <30km or >200km
#out=final_events_P.lst
out=final_events_P_XJ1997.lst
#out=final_events_P_XE2005-2007.lst
#out=final_events_S_TO2013-2015.lst

# Central of area used to calculate baz of events
if false;then
###@@@@ first time @@@@@###
#*********** sort the list by baz **********
cat /dev/null > temp
#for cl in d0-30 d30-200 d200-10000;do
#for cl in P_XJ1997_d0-1000;do
for cl in P_XJ1997;do
#for cl in P_XE2005-2007_d0-1000;do
#for cl in S_TO2013-2015;do
  IFS=$'\n'
  for line in `cat final_events_${cl}.lst`;do
    evtpath=`echo $line |awk '{print $1}'` 
    pspath=`echo $evtpath |awk -F/ '{printf"%s/%s",$1,$2}'`
    evnm=`echo $evtpath |awk -F/ '{print $3}'`
    lon=`echo $line |awk '{print $2}'` 
    lat=`echo $line |awk '{print $3}'` 
    dep=`echo $line |awk '{print $4}'` 
    dist=`~/progs/sactools_c/distaz 36.5 -119.5 $lat $lon |tail -1 |awk '{print $1}' `
    baz=`~/progs/sactools_c/distaz 36.5 -119.5 $lat $lon |tail -2 |awk '{print $2}' `
    ~/progs/TauP-2.4.5/bin/taup_time -mod prem -h $dep  -ph P -deg $dist >taup.out
    inc_ang=`cat taup.out |sed -n '6p' |awk '{print $7}'`  ## Incident angle
 
    echo $evtpath $lon $lat $dep $baz $inc_ang >> temp
  done
done
#sort -k 5 -n temp >final_events_P_XJ1997_d0-1000.lst
sort -k 5 -n temp >final_events_P_XJ1997_d0-1000.lst
#sort -k 5 -n temp >final_events_P_XE2005-2007_d0-1000.lst
#sort -k 1 -n temp >final_events_S_TO2013-2015_d0-1000.lst
exit
#********************************************
i=1
cat /dev/null >$out
#for cl in d0-10000;do
#for cl in P_XJ1997_d0-1000;do
 for cl in P_XE2005-2007_d0-1000;do
  IFS=$'\n'
  for line in `cat final_events_${cl}.lst`;do
    evtpath=`echo $line |awk '{print $1}'` 
    pspath=`echo $evtpath |awk -F/ '{printf"%s/%s",$1,$2}'`
    evnm=`echo $evtpath |awk -F/ '{print $3}'`
    lon=`echo $line |awk '{print $2}'` 
    lat=`echo $line |awk '{print $3}'` 
    dep=`echo $line |awk '{print $4}'` 
    baz=`echo $line |awk '{print $5}'` 
    inc_ang=`echo $line |awk '{print $6}'` 

    echo "===> baz,ind = $baz,$inc_ang" 

    gv -scale=0.6  $pspath/seismo_${evnm}.ps &
    echo "Save this event or Not? Yes(y)No(n)"
    read issave
    if [ $issave == "y" ];then
       echo $evtpath $lon $lat $dep  $baz $inc_ang >>$out
    fi
    let i=i+1
    if [ $i -eq 50 ];then
       killall gv
    fi

  done
done
mv $out ${out}.pre
killall gv
fi
###@@@@ Second time to check and modify @@@@@###
if true;then
cat /dev/null >$out
IFS=$'\n'
for line in `cat $out.pre`;do
    evtpath=`echo $line |awk '{print $1}'`
    pspath=`echo $evtpath |awk -F/ '{printf"%s/%s",$1,$2}'`
    evnm=`echo $evtpath |awk -F/ '{print $3}'`
    lon=`echo $line |awk '{print $2}'`
    lat=`echo $line |awk '{print $3}'`
    dep=`echo $line |awk '{print $4}'`
    baz=`echo $line |awk '{print $5}'`
    inc_ang=`echo $line |awk '{print $6}'`

    echo "===>dep,baz,ind = $dep km $baz, $inc_ang" 

    gv -scale=0.6  $pspath/seismo_${evnm}.ps &
    echo "Save this event or Not? Yes(y)No(n)"
    read issave
    if [ $issave == "y" ];then
       echo $evtpath $lon $lat $dep $baz $inc_ang >>$out
    fi
    let i=i+1
    if [ $i -eq 10 ];then
       killall gv
    fi
done
killall gv
fi
