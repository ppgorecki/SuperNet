#!/usr/bin/fish

 
set --local options 'h/help' 't/contracttest' 'r/rtcount=!_validate_int --min 0' 'l/lfcount=!_validate_int --min 2' 'L/retusageleft=!_validate_int --min 0' 'R/retusageright=!_validate_int --min 0' 'd/draw' 'n/net=' 'g/genetree=' 'D/dc' 'P/dptest' 'B/bbtest=' 'H/hctest='

argparse $options -- $argv; or set _flag_help 1

if set --query _flag_help
    printf "Usage: $0 [OPTIONS]\n\n"
    printf "Options:\n"
    printf "  -h/--help     Prints help and exits\n"
    printf "  -t/--contracttest       Prints help and exits\n"
    printf "  -r/--rtcount=NUM  Ret. count (minimum 0, default 5)\n"
    printf "  -l/--lfcount=NUM  Leaf count (minimum 2, default 2)\n"
    printf "  -d/--draw=NUM  Draw\n"
    printf "  -L/--retusage=NUM  Ret usage left\n"
    printf "  -R/--retusage=NUM  Ret usage right\n"
    printf "  -g/--genetree=TREE  Gene tree\n"    
    printf "  -D/--dc=TREE  DC test\n"    
    printf "  -P/--dptest  DP/contr test\n" 
    printf "  -B/--bbtest  BB test\n" 
    printf "  -H/--hctest  HC test\n"       
    exit 0
end

set --query _flag_rtcount; or set --local _flag_rtcount 5
set --query _flag_lfcount; or set --local _flag_lfcount 2
set --query _flag_retusageleft; or set --local _flag_retusageleft 0
set --query _flag_retusageright; or set --local _flag_retusageright 0

if [ $_flag_contracttest ];
    
    set pw (echo "print(2 ** $_flag_rtcount - 1)" | python3 )     
    set cnt 1 
    while true
        set x (random 0 $pw)
        set y (echo "print($x ^ $pw)" | python3 )
        echo "./supnet -R$_flag_rtcount -A$_flag_lfcount -r1 -eX1 0 0"
        set N (./supnet -R$_flag_rtcount -A$_flag_lfcount -r1 -eX1 0 0)

        echo $cnt './supnet -n "'$N'"' -eX $x $y
        if ! ./supnet -n $N -eX $x $y > x.txt; break; end
        tail -1 x.txt             
        set cnt (math $cnt + 1)
    end
    testcontr.fish -N (head -1 x.txt) -L $x -R $y 
    exit 0
end

if [ $_flag_draw ];    
    set R cx
    #echo "./supnet -R$_flag_rtcount -A$_flag_lfcount -r1 -e1n"    
    set --query _flag_net; or set --local _flag_net (./supnet -R$_flag_rtcount -A$_flag_lfcount -r1 -e1n)        
    ./supnet -n $_flag_net -e1Xng $_flag_retusageleft $_flag_retusageright     
    echo "./supnet $_flag_genetree -n $_flag_net -e1Xng $_flag_retusageleft $_flag_retusageright "
    mv contr.dot $R.dot 
    dot -Tpdf $R.dot -o $R.pdf
end

if [ $_flag_dc ];    
    set R cx    
    #echo $_flag_rtcount
    set --query _flag_net; or set --local _flag_net (./supnet -R$_flag_rtcount -A$_flag_lfcount -r1 -e1n)        
    set --query _flag_genetree; or echo gene tree expected && exit -2
    
    ./supnet -g $_flag_genetree -n $_flag_net -eXng $_flag_retusageleft $_flag_retusageright  | tail -2 > x.log
    set retmindc ( cat x.log | head -1 | cut -f2 -d: )
    set contrnet ( cat x.log | tail -1 )
    set retmindcc (./supnet -g $_flag_genetree -n $contrnet -ed)    
    echo $retmindc $retmindcc $contrnet $_flag_net
    if [ $retmindcc != $retmindc ];         
        echo
        echo "==============="
        echo 

        echo "./supnet -g \"$_flag_genetree\" -n \"$_flag_net\" -eXng $_flag_retusageleft $_flag_retusageright"
        ./supnet -g $_flag_genetree -n $_flag_net -eXng $_flag_retusageleft $_flag_retusageright

        echo "==============="

        echo "./supnet -g \"$_flag_genetree\" -n \"$contrnet\" -ed"
        ./supnet -g $_flag_genetree -n $contrnet -ed

        echo "===============" embnet        
        echo "embretnet/embnet.py -g $_flag_genetree -n $contrnet -pe"
        embretnet/embnet.py -g $_flag_genetree -n $contrnet -pe
        
        exit -2; 
    end


end

if [ $_flag_dptest ];
    
    set pw (echo "print(2 ** $_flag_rtcount - 1)" | python3 )     
    set cnt 1 
    while true
        set x (random 0 $pw)
        set y (random 0 $pw)
        set y (echo "print($y ^($y & $x))" | python3 )
        # set y (echo "print($x ^ $pw)" | python3 )

        set N (./supnet -R$_flag_rtcount -A$_flag_lfcount -r1 -eX1 0 0)

        echo TEST $cnt  "./supnet -eX $x $y -n \"$N\""
        if ! testcontr.fish -D -L$_flag_retusageleft -R$_flag_retusageright -l$_flag_lfcount -g $_flag_genetree; 
            break; 
        end
        
        #if ! ./supnet -n $N -eX $x $y > x.txt; break; end
        #tail -1 x.txt             
        set cnt (math $cnt + 1)
    end
    #testcontr.fish -N (head -1 x.txt) -L $x -R $y 
    exit 0
end

if [ "$_flag_bbtest" = 'r' ];
    set cnt 0
    while testcontr.fish -B x; 
        echo OK $cnt
        set cnt (math $cnt + 1) 
    end
    exit 0
end

# supnet only
if [ "$_flag_bbtest" = 's' ];
    set net (supnet -r1 -R10 -A20 -e1n)
    set gtree (supnet -r1 -R0 -A20 -en)
    set cnt 0
    while supnet -n $net -g $gtree -eY | tail -1
        echo OK $cnt
        set cnt (math $cnt + 1) 
    end
    exit 0
end

# supnet only
if [ "$_flag_bbtest" = 'c' ];
    set MX 200
    set net (supnet -r$MX -R5 -A20 -en)
    set gtree (supnet -r$MX -R0 -A20 -en)
    
    time for i in (seq 1 $MX);
        supnet -n $net[$i] -g $gtree[$i] -eY | tail -1        
    end > tmp/supnet.log
    time for i in (seq 1 $MX);
        embretnet/embnet.py -n $net[$i] -g $gtree[$i] -pb | tail -1        
    end > tmp/embnet.log

    diff tmp/supnet.log tmp/embnet.log
    exit 0
end       



if [ $_flag_bbtest ];    

     if [ $_flag_bbtest = l ];
        set net (cat tmp/bbnet.txt)
        set gtree (cat tmp/bbgtree.txt)
     else
        set net (supnet -r1 -R10 -A10 -e1n)
        set gtree (supnet -r1 -R0 -A10 -en)
        #set gtree (embretnet/embnet.py -g rand:10:0 -pg)
        #set net (embretnet/embnet.py -n rand:13:0  -pn)
        #set gtree (embretnet/embnet.py -g rand:13:0 -pg)
        echo $net > tmp/bbnet.txt
        echo $gtree > tmp/bbgtree.txt
     end

     echo $net $gtree

     embretnet/embnet.py -n $net -g $gtree -pb > tmp/e.log
     cat e.log

     set ec (tail -1 tmp/e.log)

     echo TEST 
     
     #supnet -n $net -g $gtree -eY 

     set es (supnet -n $net -g $gtree -eY | tail -1 )

     echo "supnet -n \"$net\" -g \"$gtree\" -eY"

     echo $es $ec

     test "$es" = "$ec" ; or exit -1

     exit 0
end



if [ $_flag_hctest ];    

     if [ $_flag_hctest = l ];
        set net (cat tmp/hcnet.txt)
        set gtree (cat tmp/hcgtree.txt)
     else
        set net (supnet -r1 -R20 -A10 -e1n)
        supnet -r3 -R0 -A10 -en > tmp/hcgtree.txt
        echo $net > tmp/hcnet.txt        
     end

     echo $net $gtree

     echo "DP&BB"
     supnet -n $net -G tmp/hcgtree.txt -oT1t -CDC
     mv odt.dat tmp/dpbb_odt.dat

     echo "Exhaustive"
     supnet -n $net -G tmp/hcgtree.txt -oT1te -CDC
     mv odt.dat tmp/naive_odt.dat

     

     cat tmp/naive_odt.dat
     cat tmp/dpbb_odt.dat
     

     exit 0
end
