#!/bin/bash
for i in 1.15 1.25 1.35 1.45 1.55 1.65 1.75 1.85 1.95 2.05
do
 make T=$i Data$i.dat &
done
