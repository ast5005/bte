#BTEmain.jl

module BTEmain
using PyPlot

using constants
dos=Float64{}
w=collect(0:0.1:1.0)
include("DOSfunctions.jl")
Emin=0.0
Einterval=0.01
Emax=1.0
#using types
function simulate()
  #x=testFunction(12354635)
  #println("Result=",x)
    
    parBand0=parBand(1.0*me,0.0,0.0)
    dos=getDOS_E(parBand0,Emin,Emax,Einterval)
    println("DOS matrix:\n",dos,"\n")
    E=collect(Emin:Einterval:Emax)
    println("Energy Levels:\n",E,"\n")
    plot(x=E,y=dos)
end
simulate()

end