# Formation of wide exoKuiper belts from migrating dust traps
Python wrapper for DustPy v2 to simulate evolution of planetesimal formation resulting in wide debris discs.

Authors: Elle Miller & Sebastian Marino (2021)

## File Tree

There are lots of files I created over the project, but you will most likely only end up using a few. The bare minimum files to run a simulation are 
- main.py
- functions.py
- functionsPlanFormation.py
- functionsMovingBump.py
- bumpParams.py

The most important file here to familiarise yourself with is **main.py**, which is the 'control' file for running the program and sets all the initial conditions. 

To run a simulation:

``$ python main.py -flag1 val1 -flag2 val2`` etc.

For example, to run a simulation of the prime case (alpha=1e-3, A=10, v=100%, pos=90au) and store the results in directory number 3 you would do:

``$ python main.py -z 3 -n 31 -a 1e-3 -b 10 -v 1 -p 90``

- The _plot_ folder contains scripts to plot data. Each script plots a different thing, and so you need to check the file and its arguments (at the bottom of the file) to see what the plotting options are. For example, ``python plotSDR.py -z 1`` plots the SD profiles of the last snapshot in directory '1'. There's loads of functionality here for you to use, such as creating movies (``plotSim.py -s 1``), but many instructions may not be clearly documented and so msg me if you get stuck.

- The _run_ folder contained my bash files to run a suite of jobs at once in Slurm. If your computational cluster also uses Slurm then I can show you how to use these, but most likely it will be something different. I will note here that in a project that involves running many simulations, you will be far better off in the long run going to the effort of automating as much you can from the start.

In cleaning up this repo (Feb 2022) I moved all the plotting files into the **plot** folder and run scripts into the **run** folder, without updating the necessary directory code in the files. You will need to go through and add all of your own directories anyway, so remember to also to this e.g. ../main.py
