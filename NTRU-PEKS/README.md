# PEKS

PEKS for protected soft-biometric keyword search
Installation & build guide, updated 11/04/2021

Setup:
- Check if you have GMP installed in /usr/local/libntl
- If not, install GMP version 6.2.0 (not the newest version 6.2.1!) following this guide https://libntl.org/doc/tour-gmp.html
- leave out the --prefix, just type ./configure, this will install GMP in /usr/local
- in the last step, type sudo make install instead of just make install (otherwise you will get an error)
- Install NTL library as described in https://libntl.org/doc/tour-unix.html (it will also be installed in /usr/local)
- make check takes quite a while, time to grab a coffee
- Install CodeLite IDE
- Clone this repo (or https://github.com/Rbehnia/NTRUPEKS for original) into ~/$USER/Documents/MyCode
- Rename directory "peks" (or "NTRUPEKS" for orginal) to "NTRU-PEKS"
- Open CodeLite, create "new workspace" named "NTRU-PEKS" under ~/$USER/Documents/MyCode/NTRU-PEKS and uncheck the "Create the workspace under a seperate directory" box (workspace path should be ~/$USER/Documents/MyCode/NTRU-PEKS/NTRU-PEKS.workspace)
- After this, the NTRU-PEKS.workspace file should be in the same directory as the .cc files, the makefile and the project file
- Choose "Add existing project" to workspace, choose NTRU-PEKS.project file in ~/$USER/Documents/MyCode/NTRU-PEKS
If "add existing project" does not work:
- Create new project in the workspace (call it "NTRU-PEKS"), don't forget to add virtual directory "src", and add all .cc and .h files manually
Finally:
- right-click the project, go to settings, Compiler
- paste "-O2;-std=c++11;-Wall" into C++ Compiler Options and "-O2;-Wall" into C Compiler Options
- click "Apply"
- go to Linker
- past "-lntl -lgmp" into Linker Options
- click "Apply" and "Ok"
- set paths in Scheme.cc function "Keyword_Database_Setup" to "/home/$USER/Documents/MyCode/NTRU-PEKS/frgc_gender_10000", same for ethnicity, skintype, agegroup
- set nb_sub in params.h to 500
Now you are ready to build the project (build in Debug mode first)
- At first build, map suggested Cross GCC compiler to existing GCC or CLANG (CodeLite will offer you this, if not, you are fine)
- Enjoy
