#!/bin/bash
for i in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
do
 make T=$i Data$i.dat &
done
