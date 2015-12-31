cd statlib
f2py -m  as159 asa159.f90 -c  
f2py -m  fexact FEXACT.F90 -c
cd ..