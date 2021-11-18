#!/bin/bash
for i in 3.1 3.2 3.3 3.4 3.5 3.6 3.7 3.8 3.9 4.0
do
 make T=$i Data$i.dat &
done
