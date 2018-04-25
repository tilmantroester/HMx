MODULE constants

  IMPLICIT NONE

  !Mathematical constants
  REAL, PARAMETER :: pi=3.14159265359 !pi
  REAL, PARAMETER :: em=0.5772156649 !Euler–Mascheroni
  REAL, PARAMETER :: zero=0. !zero
  REAL, PARAMETER :: one=1. !one

  !Physical constants
  REAL, PARAMETER :: kb=1.38065e-23 !Boltzmann constant in SI  
  REAL, PARAMETER :: mp=1.6726219e-27 !Proton mass in kg
  REAL, PARAMETER :: bigG=6.67408e-11 !Gravitational constant in kg^-1 m^3 s^-2
  REAL, PARAMETER :: eV=1.60218e-19 !Electronvolt in Joules
  REAL, PARAMETER :: cm=0.01 !Centimetre in metres
  REAL, PARAMETER :: rad2deg=180./pi !Radians-to-degrees conversion

  !Cosmology
  REAL, PARAMETER :: mpc=3.086e22 ! m/Mpc
  REAL, PARAMETER :: yfac=8.125561e-16 !sigma_T/m_e*c^2 in SI
  REAL, PARAMETER :: conH0=2998. !(c/H0) in Mpc/h
  REAL, PARAMETER :: msun=1.989e30 ! kg/Msun
  REAL, PARAMETER :: critical_density=2.775d11 !Critical density at z=0 in (M_sun/h)/(Mpc/h)^3 (3*H0^2 / 8piG)
  REAL, PARAMETER :: dc0=(3./20.)*(12.*pi)**(2./3.) !Einstein-de Sitter linear collapse density ~1.686
  REAL, PARAMETER :: Dv0=18.*pi**2 !Einsten-de Sitter virialised collapse threshold ~178

  !Weirdly specific things
  !REAL, PARAMETER :: epn=0.875 !1/mu_e electron per nucleon (~8/7 for x=0.25 He mass frac)
  !REAL, PARAMETER :: mue=1.14 !mu_e nucleons per electron (~7/8 for x=0.25 He mass frac)
  REAL, PARAMETER :: fh=0.76 !Hydrogen mass fraction
  REAL, PARAMETER :: mue=2./(1.+fh) !Nucleons per electron (~1.136 if fh=0.76)
  REAL, PARAMETER :: pfac=(5.*fh+3.)/(2.*(fh+1.)) !Pressure factor (Hill & Pajer 2013; I do not understand; ~1.932 if fh=0.76)  
  
END MODULE constants
