#gfortran -pg -O3 -march=native -fbounds-check -fbacktrace rd_2d_spectral_mod.f90 test_A.f90 -I/usr/local/include -L/usr/local/lib -lfftw3 -lm -o rd_2d.exe
gfortran -pg -O3 -march=native -fbounds-check -fbacktrace rd_2d_spectral_mod.f90 test_B.f90 -I/usr/local/include -L/usr/local/lib -lfftw3 -lm -o rd_2d.exe
