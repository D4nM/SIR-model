# SIR-model
Simple graphical simulation of the diffusion of a pathogen in a confined and finite population, based on the basic S.I.R. ("Susceptible, Infected, Recovered") epidemiological model.  
Made in July 2020 for my Programming exam in university.

### _Nota Bene_
- This program was created and tested for an Ubuntu-derived Linux operating system. Other OSs have not been tested.
- This program was my very first "bigger" project, written all on my own and in a bit of a rush to meet the exam's deadline, thus the code might result messy and unorganised, especially since it's all in a single file.
  
# Running the program (Ubuntu)
1. Install SFML by following the steps present on [SFML's website](https://www.sfml-dev.org/tutorials/2.5/start-linux.php), or execute `$ sudo apt-get install libsfml-dev` directly from the terminal.
2. Download the archive and extract it in your desired location.
3. Navigate with your terminal to the directory containing the `.cpp` file.
4. Build & run:
   1. Build the program by executing in order:  
      `$ g++ SIR_model_v2.cpp -c`  
      `$ g++ SIR_model_v2.o -lsfml-graphics -lsfml-window -lsfml-system -o Build-SIR/SIR_model_v2`
   2. Run the program by executing: `$ ./Build-SIR/SIR_model_v2`
      - N.B.: Make sure to have the 'Fonts' directory in the same directory of the executable, otherwise no text will be seen during run time.  
 
Alternatively, you can build and run the program in one sweep by executing the script [`compile-link-run_v2.sh`](compile-link-run_v2.sh) from terminal (make it executable via `$ chmod u+x [scriptname].sh`).





