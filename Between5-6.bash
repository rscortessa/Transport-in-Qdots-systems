#!/bin/bash
for i in 5.1 5.2 5.3 5.4 5.5 5.6 5.7 5.8 5.9 6.0
do
 make T=$i Data$i.dat &
done
