radconv: radconv/constants.f90 radconv/simple_sat_vapor_pres.f90 radconv/dry_convection.f90 radconv/qe_moist_convection.f90 radconv/radiation.f90 radconv/radconv.f90
	gfortran -c radconv/constants.f90
	gfortran -c radconv/radiation.f90
	gfortran -c radconv/simple_sat_vapor_pres.f90
	gfortran -c radconv/dry_convection.f90
	gfortran -c radconv/qe_moist_convection.f90
	gfortran radconv/radconv.f90 constants.o radiation.o simple_sat_vapor_pres.o dry_convection.o qe_moist_convection.o -o radconv/radconv.out
	rm constants.o
	rm constants_mod.mod
	rm radiation.o
	rm radiation_mod.mod
	rm simple_sat_vapor_pres.o
	rm simple_sat_vapor_pres_mod.mod
	rm dry_convection.o
	rm dry_convection_mod.mod
	rm qe_moist_convection.o
	rm qe_moist_convection_mod.mod
	rm radconv_mod.mod
	