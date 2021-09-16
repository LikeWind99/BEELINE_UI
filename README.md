# BEELINE_UI
```
######################################################################
#N...........N.......................................................#
#....N...........N...................................................#
#.......N...........N................................................#
#..........N..........N..............................................#
#............N..........N............................................#
#..............N..........N...................NNN....................#
#...............N..........N................NN..NN...................#
#............NNNNNN..NNNNNNNNN............NN....NN...................#
#.........NNNNN....NNNNNN...NNN.NNN.....NN......NN.........NNNNN.....#
#........NNN.....NNNNNNNNNNNNNN...NNN..NN........N....NNNNN....NNNNNN#
#........N..N...NNNNNNNNNNNNNNN...NNN.N........NNNNNNN............NN.#
#.....NNNN......NNNNNNNNNNNNNN......NNN......NNN.................NNNN#
#...NNNN..........NNNNNNNNNN......N.NNN...NNN..................NNNNN.#
#..NNN.......................NNNNN.N.NN.NN................NNNNNN.....#
#.NNN...................NNNNN.N.N.N.NN.N................NN...........#
#.NN...............NNNNNNNNNNNNNNN.NN.....................NNNNN......#
#..NN.........NNNNN.................NNNN..NNNNN..............NN......#
#....NNN.NNNNN..................NNNNN...NN...NNNNNN.........NN.......#
#.......NN................NNNNNNNNN.........NN..N..NN....NNN.........#
#.......N...NNNNNNNNNNNNNNNNNNNNNNN........NNNNNNN...NNNNN...........#
#.....NN......NNNNNNNNNNNNNNNNNNNNNNN.....NN......N...N..............#
#....NN.NNNN.NNNNNNNNNNNNNNNNNNNNNNNNNNNNN.N......NNNNN..............#
#...N.NNNNNNN.......NNNNNNNNNNNNNNNNNNNNN..NN....NNNNNNN.............#
#..N....NNNNNN.........NNNNNNNNNNNNNNN.....NN......NNN...............#
#...NN.NNNNNN...........NNNNNNNNNNNN.......NN........................#
#.....NNNN..............N...............NNNNN........................#
#......................NNN...........NNNNNNNN........................#
#.....................NNNNNNNNNNNNNNNNNNNNNN.........................#
#.....................N.NNNNNNNNNNNNNNNNNNN..........................#
#.....................N.................NN...........................#
#.....................NN...............NN............................#
#......................NNNNNNNNNNNNNNNN..............................#
#.......................NNNNNNNNNNNNN................................#
#........................NNNNNNNNN...................................#
#.........................NNNN.......................................#
#.........................NNN........................................#
#........................NN..........................................#
#.......................N............................................#
#Created by shenAo     Time: 2021-08-10     Email: shen_a_o@qq.com   #
######################################################################
```

BEELINE: Benchmarking gEnE reguLatory network Inference from siNgle-cEll transcriptomic data.

BEELINE(https://github.com/Murali-group/Beeline) is a pipeline which provides tools to evaluate the performance of algorithms for the reconstruction of gene regulatory networks (GRNs) from single-cell RNAseq data. However, for "green hand", they often encounter many difficulties when using it for data analysis, because the process of BEELINE is a little complex that requires users to equip with the competence of programming. Therefore, in this project, we used python to reconstruct part of the code of BEELINE (as the time is pressing, not all have been converted yet), added some new functions, and built a user-friendly graphical interface, which greatly simplifies the single-cell data analysis process and improves production efficiency.

In the following pages, I will show you how to quickly deploy this awesome project!

Quick setup:
1. Download and install python(3.5+)(https://www.python.org/downloads/windows/) from the official website.
2. Download this Repository.
3. Use pip to install the required libraries(e.g. numpy, pandas...).
       You can use the following command in the current directory to quickly install: 
           ```pip install -r requirements.txt
           ```
4. Run the "window.py", and then the graphical interface will open.

If you are used to using anaconda, maybe the following process can help you:
1. Create a python environment to run the project
  - 1.1 Create a python virtual environment
         ```conda create -n your_project_name python=3.5
         ```
  - 1.2 Confirm whether the environment is created successfully
         ```conda info -e
         ```
  - 1.3 Activate the environment
         ```activate your_project_name
         ```
2. Change workspace to your current directory.
      For example: if you put the folder on the desktop, then execute 
         ```cd C:\Users\your_accout_name\Desktop\BEELINE
         ```
3. Install dependent packages.
  ```conda install -r requirements.txt
  ```
4. Run the "window.py", and then the graphical interface will open.

This will be a long-term project, we will continue to improve the entire project, if we have time in the future.

Finally, thank a lot for my teammates: Liu Yun, Hu Xinyun, Sun Zichun, Ren Zirui for refactoring part of the code of BEELINE, and their great help to the entire project. And last, especially thanks for the instructions and support our dear instructor Shi Qianqian gave to us, we couldn't achieve this without her help.

At the end, if you have any question, please can contact us by email.

