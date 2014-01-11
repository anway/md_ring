#!/bin/bash


touch plot3d.py

echo "from matplotlib import pyplot" > plot3d.py
echo "import pylab" >> plot3d.py
echo -e "from mpl_toolkits.mplot3d import Axes3D\n\n" >> plot3d.py
echo -e "data = [" >> plot3d.py

for ((t = 0; t <= 50; t += 5))
do
    for ((z = 10; z <= 50; z += 5))
    do
        tp=$(echo "scale=1; $t/10" | bc )
        zp=$(echo "scale=1; $z/10" | bc )
#echo $tp
#echo $zp
        echo -e "30.0\n1.0 1.0\n10.0 1.0 400.0 $zp\n12\n0.001 20.0 $tp 20.0" > in.txt
        ./md_ring_tz.out
        if [[ $t -eq "50" && $z -eq "50" ]]
        then
            echo -e "]\n" >> plot3d.py
        else
            echo ", " >> plot3d.py
        fi
    done
done
echo "x = []" >> plot3d.py
echo "y = []" >> plot3d.py
echo -e "z = []\n" >> plot3d.py
echo "for list in data:" >> plot3d.py
echo -e "\tx.append(list[0])\n\ty.append(list[1])\n\tz.append(list[2])\n" >> plot3d.py
echo "fig = pylab.figure()" >> plot3d.py
echo "ax = Axes3D(fig)" >> plot3d.py
echo "ax.scatter(x,y,z)" >> plot3d.py
echo "pyplot.show()" >> plot3d.py
