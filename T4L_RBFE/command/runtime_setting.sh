#!/bin/bash
#everyλ run for 0.5 ns, but you also need to change some reportly parameter for avoiding unnecessary error.
sed -i  's/250000/25000/'  md_*.namd # OUTPUT FREQUENCIES
sed -i  's/500000/50000/g'  md_*.namd  # alchEquilSteps
sed -i  's/500000/250000/g'  md_*.namd  # run
