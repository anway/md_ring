#!/bin/bash


touch plot.py

echo "import matplotlib" > plot.py
echo -e "import pylab\n\n" >> plot.py
echo -e "data = [" >> plot.py

for ((i = 0; i <= 1000; i += 10))
do
    echo -e "30.0\n1.0 1.0\n10.0 1.0 $i.0 1.0\n12\n0.001 20.0 0.0000000000000000000006 20.0" > in.txt
    ./md_ring_k0.out
    if [ $i -lt "1000" ]
    then
        echo ", " >> plot.py
    else
        echo -e "]\n" >> plot.py
    fi
done
echo "x = []" >> plot.py
echo -e "y = []\n" >> plot.py
echo "for list in data:" >> plot.py
echo -e "\tx.append(list[0])\n\ty.append(list[1])\n" >> plot.py
echo "matplotlib.pyplot.scatter(x,y)" >> plot.py
echo "matplotlib.pyplot.show()" >> plot.py
