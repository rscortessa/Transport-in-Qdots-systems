#!/bin/bash
for i in 2.15 2.25 2.35 2.45 2.55 2.65 2.75 2.85 2.95 3.05
do
 make T=$i Data$i.dat &
done
