#!/bin/bash
for i in 3.15 3.25 3.35 3.45 3.55 3.65 3.75 3.85 3.95 4.05
do
 make T=$i Data$i.dat &
done
