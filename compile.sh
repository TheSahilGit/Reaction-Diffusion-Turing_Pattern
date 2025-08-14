#gfortran -pg -O3 -march=native -fbounds-check -fbacktrace rd_2d_mod.f90 main_prog.f90 -o rd_2d.exe
#
#gfortran rd_2d_mod.f90 rd_2d_spectral_mod.f90 main_prog.f90 -I/usr/local/include -L/usr/local/lib -lfftw3 -lm
#

gfortran -pg -O3 -march=native -fbounds-check -fbacktrace \
    rd_2d_mod.f90 rd_2d_spectral_mod.f90 main_prog.f90 \
    -I/usr/local/include \
    -L/usr/local/lib -lfftw3 -lm \
    -o rd_2d.exe
