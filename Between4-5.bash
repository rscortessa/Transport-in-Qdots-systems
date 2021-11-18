#!/bin/bash
for i in 4.1 4.2 4.3 4.4 4.5 4.6 4.7 4.8 4.9 5.0
do
 make T=$i Data$i.dat &
done
