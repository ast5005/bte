# constants.jl
#__precompile__()
module constants

export h, hbar, q, kB, me, kBe, eps0, Nav

h=6.62607004e-34 # Planck constant in m2kg/s
hbar=1.0545718e-34 #Planck constant in m2kg/s
q=1.60217662e-19 #Coulombs charge of electron
kB=1.38064852e-23 #Boltzmann Constant m2 kg s-2 K-1
kBe=8.6173324e-5; #eV K-1
me=9.10938356e-31# kilograms
eps0=8.85e-12
Nav=6.0221415e23 # Avagadro's number
end