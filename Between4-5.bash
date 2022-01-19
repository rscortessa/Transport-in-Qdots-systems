#!/bin/bash
for i in 4.15 4.25 4.35 4.45 4.55 4.65 4.75 4.85 4.95 
do
 make T=$i Data$i.dat &
done
