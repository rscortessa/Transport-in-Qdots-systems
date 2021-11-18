#!/bin/bash
for i in 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8 2.9 3.0
do
 make T=$i Data$i.dat &
done
