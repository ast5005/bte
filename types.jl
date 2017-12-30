#  types.jl

module types
import constants
import materialconstants
export parBand, BandStruc, tau_electron_base, tau_electron_AC, tau_electron_B, parBandTx, BandStrucTx, tau_phonon_B, tau_phonon_Base

    abstract type tau_electron_Base end
#should have a value
    abstract type tau_phonon_Base end

type parBand{T}
  effMass::T
  alpha::T
  offset::T
  degen=1.0
    
end
type parBandTx
    effMass::Float64
    alpha::Float64
    offset::Float64
    funsofoffset::Vector{Function} 
    funsofeffMass::Vector{Function}
    var::Vector{Float64}
end

type BandStrucTx
    bands::Vector{parBandTx}
    var::Vector{Float64}
end

type BandStruc
    bands::Vector{parBand}
end

type tau_electron_AC <: tau_electron_Base
    tauMethods::Vector{Function}
    variables::Vector{Any}      #variable[2] is T, variable[1] is meff variable[3] is E
end
type tau_electron_B <: tau_electron_Base
    tauMethods::Vector{Function}
    variables::Vector{Any}      #variable[2] is T, variable[1] is meff variable[3] is E
end
type tau_phonon_B <: tau_phonon_Base
    tauMethods::Vector{Function}
    variables::Vector{Any}      #variable[2] is T, variable[1] is meff variable[3] is E
end
function bandTxupdate(bandTx::parBandTx)
    for method in bandTx.funsofoffset
        bandTx.offset=method(bandTx.var) 
    end
    for method in bandTx.funsofeffMass
        bandTx.effMass=method(bandTx.var) 
    end
end
function updatebnstTx(bndstTx::BandStrucTx)
    for bandTx in bndstTx.bands
        bandTx.var=bndstTx.var
        bandTxupdate(bandTx)
    end
end

function bandUpdate(bndst::BandStruc,funsofoffset::Vector{Function},funsofeffMass::Vector{Function},var::Vector{Float64})
    #if length(bndst.bands)<4
    #    error("Ecpected band structure inludes 4 bands")
    #else
        for (i,band) in enumerate(bndst.bands)
            band.offset=funsofoffset[i](var)
            band.effMass=funsofeffMass[i](var)
        end
   # end
end
function calc_tau_electron_AC(variables::Vector)
    if length(variables)<3 
        error("calc_tau_electron_AC needs 3 variables m,T,E")
    end
    return pi*constants.hbar^4*materialconstants.Cl/(sqrt(2)*variables[1]^1.5*materialconstants.Da^2*constants.kB*variables[2]*(variables[3]*constants.q)^0.5)
end

function update_tau(tau_electron_Base)
    results=Array{Any}(length(tau_electron_Base.tauMethods))
    for (i,methods) in enumerate(tau_electron_Base.tauMethods)
        methods(tau_electron_Base.variables) 
    end    
end
type tempDepProp
    TempVector::Vector{Float64}
    bndst::BandStruc
    tau_Temp::tau_electron_Base
    update_Temp::Vector{Function}
    update_All::Function    
end
function update_All(tdp::tempDepProp,Temp::Vector{Float64})
    for f in  tdp.update_Temp
        f(Temp) 
    end        
end

        
        
#types module ends
end


