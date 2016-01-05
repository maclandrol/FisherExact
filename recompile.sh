cd statlib
f2py -m  asa159 asa159.f90 -c  
f2py -m  asa205 asa205.f90 -c 
f2py -m  fexact FEXACT.F90 -c
cd ..