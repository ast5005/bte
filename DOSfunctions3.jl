
# DOSfunctions.jl
#module DOSfunctions

using constants
using types
using Optim
using Roots
using QuantEcon

##############################################################################################################
##############################################################################################################
#DOS Functions
#
#
#
const hbar3=1.1728121633258218e-102
const pi2=9.869604401089358
intgrid2=200
gridx=1000
function setDOSintgrid(x)
     intgrid2=x
     return 1
end
function setDOSgrid(x)
     gridx=x
     return 1
end
function getDOS_SingleBand_E(band::parBandTx,E::Float64)
    return getDOS_SingleBand_E(band,[E])[1]
end
function getDOS_SingleBand_E(band::parBandTx,E::Array{Float64,1})
    effMass=0.0
    n=length(E)
    dos=Array{Float64}(n) 
    alpha=band.alpha/q
    #println("here $alpha")
    for i in 1:n
     if band.effMass<0.0 
            if band.offset>E[i]
                effMass=-1*band.effMass                
                #dos[i]=sqrt(2).*sqrt((effMass)^3)./hbar^3./pi^2.*sqrt((-E[i].+band.offset)*q)
                Et=(-E[i].+band.offset)*q                                
                dos[i]=sqrt(2).*sqrt((effMass)^3)./hbar^3./pi^2.*sqrt(Et+alpha*Et.^2).*(1+2alpha*Et)
            #println("dos variable in getDOS_SingleBand_E meff<0 ",dos)
             else
                dos[i]=0.0
             end
        end
        if band.effMass>0.0 
            if E[i]>=band.offset                         
                #dos[i]=sqrt(2).*(band.effMass)^1.5./hbar^3./pi^2.*sqrt((E[i].-band.offset)*q)
                Et=(E[i].-band.offset)*q
                dos[i]=sqrt(2).*sqrt((band.effMass)^3)./hbar^3./pi^2.*sqrt(Et+alpha*Et.^2).*(1+2alpha*Et)
            #println("dos variable in getDOS_SingleBand_E meff>0 ",dos)
             else
                dos[i]=0.0
             end
        end
    end
    return dos
end
#
#
#
function getDOS_SingleBand_E_Range(band,E_lowlimit::Float64,E_highlimit::Float64,interval::Float64)  
    Erange=collect(E_lowlimit:interval:E_highlimit)           
    dos=Float64[]  
     println("dos at getDOS_E():\n",dos)
    for E in Erange         
        push!(dos,getDOS_SingleBand_E(band,E))
    end
 return dos
end
#
#
#
function getDOS_MultiBand_E_Total(bndst,E)
    totalDOS=0.0
    for band in bndst.bands
        totalDOS+=getDOS_SingleBand_E(band,E)
    end
    return totalDOS
end
#
#
#
function getDOS_MultiBand_E_Range(bndst,E_lowlimit::Float64,E_highlimit::Float64,interval::Float64)
    Erange=collect(E_lowlimit:interval:E_highlimit)  
    numofbands=length(bndst.bands)
    dos=Array{Float64}(length(Erange),numofbands)
    for (i,band) in enumerate(bndst.bands)      
        dos[:,i]=getDOS_SingleBand_E_Range(band,E_lowlimit,E_highlimit,interval)
    end
    return dos
end


#DOS Functions
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#Fermi-Dirac Functions
#
#
#
function fermiStat_Temp_Ef_E(Temp,Ef,E::Float64)
   return 1.0./(1.0+exp(q*(E-Ef)./kB./Temp))    
end
#
#
function fermiStat_Temp_Ef_E(Temp,Ef,E::Array{Float64,1}) # return dimensionless
   return 1.0./(1.0+exp(q.*(E.-Ef)./kB./Temp))    
end
#
#
#
function fermiStat_Temp_Ef_E_Range(Ef::Float64,Temp::Float64,E_lowlimit::Float64,E_highlimit::Float64,interval::Float64)
    Erange=collect(E_lowlimit:interval:E_highlimit)
    #println("E coming to fermiStat():\n",Erange)
    fermifunction=Float64[]
    for E in Erange
        push!(fermifunction,1.0./(1.0+exp(q*(E-Ef)/kB/Temp)))
    end
    return fermifunction    
end
#
#
#
function fermiDerivativeTemp_Ef_E_Range(Ef::Float64,Temp::Float64,E_lowlimit::Float64,E_highlimit::Float64,interval::Float64)
    Erange=collect(E_lowlimit:interval:E_highlimit)
    fermiD=Float64[]
    for E in Erange        
        push!(fermiD,fermiDerivativeTemp_Ef_E(Ef,Temp,E))
    end
    return fermiD    
end
#
#
#
function fermiDerivativeTemp_Ef_E(Ef::Float64,Temp::Float64,E::Float64)#return 1/eV    
    return  fermiDerivativeTemp_Ef_E(Ef,Temp,[E])[1]
end
#
#
#
function fermiDerivativeTemp_Ef_E(Ef::Float64,Temp::Float64,E::Array{Float64,1})#return 1/eV
    Q=exp((E.-Ef)./kBe./Temp)
    return  -Q./kB./Temp./(1.0+Q).^2 
end
#
#
#
function Fermilevel_n2(numberofn,bndst,Temp,xmax)
    C=Float64[]
    V=Float64[]
    for band in bndstTx.bands
        if band.effMass<0
            push!(V,band.offset)
        end
         if band.effMass>0
            push!(C,band.offset)
        end
    end
    equation(x)=abs(numberofn*1e6-NumofnMultiBand(bndst,x,Temp,xmax))
    if numberofn<-1e12
        EfcalcM=optimize(equation,minimum(C)-20kBe*Temp,10.0,method=:golden_section)
    elseif numberofn>1e12
        EfcalcM=optimize(equation,0.0,maximum(V)+20kBe,method=:golden_section)
    else
        error("numberofn cannot be so small >1e12 or <-1e12")
    end
    #return (EfcalcM.minimum,EfcalcM.converged)
    return EfcalcM.minimum
end
function Fermilevel_n(numberofn::Float64,bndst::BandStrucTx,Temp,xmax)
    C=Float64[]
    V=Float64[]
    for band in bndstTx.bands
        if band.effMass<0
            push!(V,band.offset)
        end
         if band.effMass>0
            push!(C,band.offset)
        end
    end
    equation(x)=numberofn*1e6-NumofnMultiBand2(bndst,x,Temp,xmax)
    if numberofn<-1e12
        #println("numberofn<-1e12")
        EfcalcM=fzero(equation,minimum(C)-20kBe*Temp,10.0)#
    elseif numberofn>1e12
        #println("numberofn<-1e12")
        EfcalcM=fzero(equation,0.0,maximum(V)+20kBe)
    else
        error("numberofn cannot be so small >1e12 or <-1e12")
    end
    #return (EfcalcM.minimum,EfcalcM.converged)
    return EfcalcM
end
#Fermi-Dirac Functions
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#Velocity
function velocity(band,E_lowlimit::Float64,E_highlimit::Float64,interval::Float64)
   Erange=collect(E_lowlimit:interval:E_highlimit)
    #println("E coming to fermiStat():\n",Erange)
    vel=Float64[]
    for E in Erange
        push!(vel,sqrt(square_velocity_E(band,E)))
    end
    return vel    
end
function velocity_E(band,E::Float64)       
    return vel= E>0 ? sqrt(2./3./abs(band.effMass).*E.*(1+band.alpha.*E)./(1+2*band.alpha.*E)^2): -sqrt(2./3./abs(band.effMass).*-E.*(1+band.alpha.*-E)./(1+2*band.alpha.*-E)^2)
end
function square_velocity_E(band,E::Float64) 
    vel=0.0
    if band.effMass>=0        
        vel= 2./3./band.effMass.*(E-band.offset).*q.*(1+band.alpha.*E)./(1+2*band.alpha.*E)^2     
    elseif band.effMass<0        
        vel= -2./3./band.effMass.*(-E+band.offset).*q.*(1+band.alpha.*E)./(1+2*band.alpha.*E)^2
    end
return vel   
end
function square_velocity_E2(band,E::Float64) 
    vel=0.0
    if band.effMass>0.0 
        if E>=band.offset   
            Et=(E-band.offset)
            vel= 2./3./band.effMass.*Et*q.*(1+band.alpha.*Et)./(1+2*band.alpha.*Et)^2          
        else
            vel=0.0
        end
    elseif band.effMass<0  
        if band.offset>E                          
            #dos[i]=sqrt(2).*sqrt((effMass)^3)./hbar^3./pi^2.*sqrt((-E[i].+band.offset)*q)
            Et=(-E.+band.offset)
            vel= -2./3./band.effMass.*Et.*q.*(1+band.alpha.*Et)./(1+2*band.alpha.*Et)^2
        else
            vel=0.0
        end
    end
return vel   
end
function square_velocity_E(band,E::Array{Float64})
    vel=Array{Float64}(length(E))
    for (i,Ex) in enumerate(E)
        vel[i]=square_velocity_E(band,Ex)    
    end
    return vel
    
end
#
#
#
function square_velocity_E_Range(band,E_lowlimit::Float64,E_highlimit::Float64,interval::Float64)
   Erange=collect(E_lowlimit:interval:E_highlimit)
    #println("E coming to fermiStat():\n",Erange)
    vel=Float64[]
    for E in Erange
        push!(vel,square_velocity_E(band,E))
    end
    return vel    
end
#Velocity
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#Number of Electrons
function Numofn(band,Ef,Temp::Float64,xmax::Float64)
    if band.effMass>=0.0         
        integrandp_Numofn(E)=q*getDOS_SingleBand_E(band,E).*(-1*fermiStat_Temp_Ef_E(Temp,Ef,E))
        #a=quadgk(integrand,band.offset,band.offset+20kBe*Temp)
        #println("calculating Numofn band.effMass>0")
        a=quadgk(integrandp_Numofn,0,100.0)        
    else
        min=band.offset-20kBe*Temp<0 ? 0.0 : band.offset-20kBe*Temp 
        integrandn_Numofn(E)=q*getDOS_SingleBand_E(band,E).*(1-fermiStat_Temp_Ef_E(Temp,Ef,E))
        #println("min is $min and band is ",band.offset)
        a=quadgk(integrandn_Numofn,min,band.offset)        
    end    
    return a[1]
end

function Numofn2(band,Ef,Temp::Float64,xmax::Float64)
    #gridx=100
    if band.effMass>=0.0         
        integrandp_Numofn2(E)=q*getDOS_SingleBand_E(band,E).*(-1*fermiStat_Temp_Ef_E(Temp,Ef,E))
        #a=quadgk(integrand,band.offset,band.offset+20kBe*Temp)
        nodes, weights = qnwlege(gridx, band.offset, band.offset+100kBe*Temp)
        #nodes, weights = qnwlege(100, Ef-50kBe*Temp, Ef+50kBe*Temp)
        #nodes, weights = qnwlege(10000, 0,100)
        a= do_quad(integrandp_Numofn2,nodes, weights)
    else
        min=band.offset-100kBe*Temp<0 ? 0.0 : band.offset-100kBe*Temp 
        #min=Ef-300kBe*Temp<0 ? 0.0 : Ef-300kBe*Temp 
        integrandn_Numofn2(E)=q*getDOS_SingleBand_E(band,E).*(1-fermiStat_Temp_Ef_E(Temp,Ef,E))
        #println("min is $min and band is ",band.offset)
        #a=quadgk(integrandn_Numofn2,min,band.offset)
        #nodes, weights = qnwlege(100, min, band.offset)
        #nodes, weights = qnwlege(10000, 0, 100)
        nodes, weights = qnwlege(gridx, min, band.offset)
        a= do_quad(integrandn_Numofn2,nodes, weights)
        
    end    
    return a[1]
end

function NumofnDerivativeEf(band,Ef,Temp::Float64,xmax::Float64)
    if band.effMass>0.0
        min=Ef-10kBe*Temp<0 ? 0.0 : Ef-10kBe*Temp 
        integrandp_NumofnDerivativeEf(E)=q*getDOS_SingleBand_E(band,E).*(fermiDerivativeTemp_Ef_E(Temp,Ef,E))
        a=quadgk(integrandp_NumofnDerivativeEf,min,Ef+10kBe*Temp)
    else
        min=Ef-10kBe*Temp<0 ? 0.0 : Ef-10kBe*Temp 
        integrandn_NumofnDerivativeEf(E)=q*getDOS_SingleBand_E(band,E).*(-1+fermiDerivativeTemp_Ef_E(Temp,Ef,E))
        a=quadgk(integrandn_NumofnDerivativeEf,min,Ef+10kBe*Temp)
    end    
    return a[1]
end
function NumofnMultiBand(bndst::BandStrucTx,Ef,Temp::Float64,xmax::Float64)
    totalnumofn=0.0
    for band in bndst.bands
        totalnumofn=totalnumofn+Numofn(band,Ef,Temp,xmax) 
    end
    return totalnumofn 
end
function NumofnMultiBand2(bndst::BandStrucTx,Ef,Temp::Float64,xmax::Float64)
    totalnumofn=0.0
    for band in bndst.bands
        totalnumofn=totalnumofn+Numofn2(band,Ef,Temp,xmax) 
    end
    return totalnumofn 
end

function NumofnMultiBandDerivativeEf(bndst,Ef,Temp::Float64,xmax::Float64)
    totalnumofn=0.0
    for band in bndst.bands
        totalnumofn=totalnumofn+Numofn(band,Ef,Temp,xmax) 
    end
    return totalnumofn 
end

#Number of Electrons
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#Relaxation Time
#
#
#
function get_tau(tau_electron_Base)
    tau=0.0
    for (i,methods) in enumerate(tau_electron_Base.tauMethods)
        tau+=methods(tau_electron_Base.variables) 
    end
    return tau
end
function get_tau(tau_electron_Base,E)
    tau=0.0
    tau_electron_Base.variables[3]=E
    for (i,methods) in enumerate(tau_electron_Base.tauMethods)
        a=methods(tau_electron_Base.variables)         
        tau+=1./methods(tau_electron_Base.variables) 
    end
    return 1./tau
end
function get_tau_total(tau_electron_Base,E)
    taud=0.0
    #tau_electron_Base.variables[1]=band.effMass
    tau_electron_Base.variables[3]=E
    for (i,methods) in enumerate(tau_electron_Base.tauMethods)
        taud+=1/methods(tau_electron_Base.variables) 
    end
    return 1/taud
end
function get_tau_phonon(tau_phonon_Base,x)
    tau=0.0
    tau_phonon_Base.variables[1]=x
    for (i,methods) in enumerate(tau_phonon_Base.tauMethods)
        a=methods(tau_phonon_Base.variables)         
        tau+=1./methods(tau_phonon_Base.variables) 
    end
    return 1./tau
end
function tau_integral(tau_electron_Base)   
    integrand(x)=q*getDOS_SingleBand_E(tau_electron_Base.variables[6],x)*-fermiDerivativeTemp_Ef_E(tau_electron_Base.variables[5],tau_electron_Base.variables[2],x)  
    minint=tau_electron_Base.variables[5]-10kBe*tau_electron_Base.variables[2]<0 ? 0.0 : tau_electron_Base.variables[5]-10kBe*tau_electron_Base.variables[2] 
    return int=quadgk(integrand,minint,Ef+10kBe*tau_electron_Base.variables[2])[1]
end
#Relaxation Time
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#Differential Conductivity
#
#
#
function sigmaD(tau_electron::types.tau_electron_Base,band,E,Ef,Temp)   
    tau_electron.variables[1]=band.effMass
    tau_electron.variables[5]=Ef
    tau_electron.variables[6]=band
    tau=get_tau(tau_electron,E)
    v=(square_velocity_E(band,E))
    temp=q*q*q*v.*tau.*getDOS_SingleBand_E(band,E).*-fermiDerivativeTemp_Ef_E(Ef,Temp,E) 
    #println(temp)
    return temp#q.*q.*v.*tau.*q.*getDOS_SingleBand_E(band,E).*-fermiDerivativeTemp_Ef_E(Ef,Temp,E)  
end
#
#
#
function sigmaD_Range(tau_electron::types.tau_electron_Base,band,Ef,Temp,E_lowlimit,E_highLimit,interval)
    Erange=collect(E_lowlimit:interval:E_highlimit)
    sigmaD=Float64[]
    for E in Erange
        push!(sigmaD,q*q*get_tau(tau_electron,E)*(square_velocity_E(band,E))*q*getDOS_SingleBand_E(band,E)*-fermiDerivativeTemp_Ef_E(Ef,Temp,E)  )
    end
    return sigmaD
end
#
#
#
function sigmaD_constanttau_Range(tau,band,Ef,Temp,E_lowlimit,E_highLimit,interval)
    Erange=collect(E_lowlimit:interval:E_highlimit)
    sigmaD=Float64[]
    for E in Erange
        push!(sigmaD,q*q*tau*(square_velocity_E(band,E))*q*getDOS_SingleBand_E(band,E)*-fermiDerivativeTemp_Ef_E(Ef,Temp,E))
        #push!(sigmaD,tau*(square_velocity_E(band,E)))
    end
    return sigmaD
end
#
#
#

#Differential Conductivity
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#Electrical Conductivity
function sigma(tau_electron::types.tau_electron_Base,band,Ef,Temp)    
    if band.effMass>0                
        integrandp_sigma(E)=q*sigmaD(tau_electron::types.tau_electron_Base,band,E,Ef,Temp) 
        nodes, weights = qnwlege(intgrid2, band.offset,band.offset+20kBe*Temp)
        a= do_quad(integrandp_sigma,nodes, weights)
        return a#quadgk(integrand,band.offset,band.offset+20kBe*Temp)[1]
    elseif band.effMass<0
        min=band.offset-20kBe*Temp<0 ? 0.0 : band.offset-20kBe*Temp    
        integrandn_sigma(E)=q*sigmaD(tau_electron::types.tau_electron_Base,band,E,Ef,Temp) 
        nodes, weights = qnwlege(intgrid2, min,band.offset)
        a= do_quad(integrandn_sigma,nodes, weights)
        return a#quadgk(integrand,min,band.offset)[1]
    else
        error("Band EffMass cannot be =0")
        return 0.0
    end
end
#
#
#
function sigma_constanttau(tau,band,Ef,Temp)
    min=Ef-10kBe*Temp<0 ? 0.0 : Ef-10kBe*Temp 
    integrandp_sigma_constanttau(E)=q*q*q*tau*(square_velocity_E(band,E))*getDOS_SingleBand_E(band,E)*-fermiDerivativeTemp_Ef_E(Ef,Temp,E)  
    #integrand(E)=(square_velocity_E(band,E))  
    return  quadgk(integrandp_sigma_constanttau,min,(Ef+10kBe*Temp))[1]
end
#
#
#
function sigma_Multiband(tau_electron::types.tau_electron_Base,bndst,Ef,Temp)
    sigmaTotal=0.0
    sigmax=0
    bandpre=bndst.bands[1]
    #sigmaTotal=sigma(tau_electron,band,Ef,Temp)
    for (i,band) in enumerate(bndst.bands)
        if i==1
            sigmax=sigma(tau_electron,band,Ef,Temp)
        else
            sigmax= band==bandpre ? sigmax : sigma(tau_electron,band,Ef,Temp) 
        end
        sigmaTotal+= sigmax
        bandpre=band
    end
    return sigmaTotal
end
function sigma_constanttau_Multiband(tau,bndst,Ef,Temp)
    sigmaTotal=0.0
    for band in bndst.bands
        sigmaTotal+= sigma_constanttau(tau,band,Ef,Temp)
    end
    return sigmaTotal
end
#Electrical Conductivity
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#Seebeck
#
#
#
function seebeck(tau_electron::types.tau_electron_Base,iband,bndst,Ef,Temp)
    seebeckTotal=0.0
    seebeckTotalDenominator=0.0
    if iband.effMass<0 
        for band in bndst.bands
            seebeckTotalDenominator+=band.effMass<0 ? sigma(tau_electron::types.tau_electron_Base,band,Ef,Temp):0.0
        end  
    end
    if iband.effMass>=0 
        for band in bndst.bands
            seebeckTotalDenominator+=band.effMass>=0 ? sigma(tau_electron::types.tau_electron_Base,band,Ef,Temp):0.0
        end  
    end   
    s=seebeck_Nominator(tau_electron::types.tau_electron_Base,iband,Ef,Temp)
    seebeckTotal= (s==0 ? 0.0:s/seebeckTotalDenominator)
     
    return seebeckTotal
end
function seebeck_sigma(tau_electron::types.tau_electron_Base,iband,bndst,Ef,Temp)
    seebeckTotal=0.0
    seebeckTotalDenominator=0.0
    if iband.effMass<0 
        for band in bndst.bands
            seebeckTotalDenominator+=band.effMass<0 ? sigma(tau_electron::types.tau_electron_Base,band,Ef,Temp):0.0
        end  
    end
    if iband.effMass>=0 
        for band in bndst.bands
            seebeckTotalDenominator+=band.effMass>=0 ? sigma(tau_electron::types.tau_electron_Base,band,Ef,Temp):0.0
        end  
    end      
    return seebeckTotalDenominator
end
#
#
#
function seebeck_constanttau(tau,iband,bndst,Ef,Temp)
    seebeckTotal=0.0
    seebeckTotalDenominator=0.0
    nominator=0.0
    if iband.effMass<0 
        for band in bndst.bands
            seebeckTotalDenominator+=band.effMass<0 ? sigma_constanttau(tau,band,Ef,Temp):0.0
        end
       # for band in bndst.bands
        nominator=seebeck_constanttau_Nominator(tau,iband,Ef,Temp)
        seebeckTotal=nominator==0 ? 0.0 : nominator/seebeckTotalDenominator
       # end 
    end
    if iband.effMass>=0 
        for band in bndst.bands
            seebeckTotalDenominator+=band.effMass>=0 ? sigma_constanttau(tau,band,Ef,Temp):0.0
        end 
      #  for band in bndst.bands
        nominator=seebeck_constanttau_Nominator(tau,iband,Ef,Temp)
        seebeckTotal=nominator==0 ? 0.0 : nominator/seebeckTotalDenominator 
      #  end 
    end       
    return seebeckTotal#seebeckTotalDenominator#
end
#
#
#
function seebeck_Nominator(tau_electron::types.tau_electron_Base,band,Ef,Temp)    
    Su=1.0
     if band.effMass>0                
        integrandseebeck_p(E)=sigmaD(tau_electron,band,E,Ef,Temp).*(E-Ef).*q
        nodes, weights = qnwlege(intgrid2,set,band.offset+20kBe*Temp)
        a= do_quad(integrandseebeck_p,nodes, weights)
        Su=a#quadgk(integrandseebeck,band.offset,band.offset+20kBe*Temp)[1]
        return -1/q/Temp*(Su)         
    elseif band.effMass<0
        min=band.offset-20kBe*Temp<0 ? 0.0 : band.offset-20kBe*Temp    
        integrandseebeck_n(E)=sigmaD(tau_electron,band,E,Ef,Temp).*(E-Ef).*q 
        nodes, weights = qnwlege(intgrid2.offset)
        a= do_quad(integrandseebeck_n,nodes, weights)
        Su=a#quadgk(integrandseebeck,min,band.offset)[1]
        return -1/q/Temp*(Su) 
    else
        error("Band EffMass cannot be =0")
        return 0.0
    end
    
    
end
#
#
#
function seebeck_constanttau_Nominator(tau,band,Ef,Temp)
    min=Ef-10kBe*Temp<0 ? 0.0 : Ef-10kBe*Temp 
    integrandseebeck(E)=q*q*q.*tau.*(square_velocity_E(band,E)).*q.*getDOS_SingleBand_E(band,E).*-fermiDerivativeTemp_Ef_E(Ef,Temp,E) *(E-Ef)
    integrandsigma(E)=q*q.*tau.*(square_velocity_E(band,E)).*q.*getDOS_SingleBand_E(band,E).*-fermiDerivativeTemp_Ef_E(Ef,Temp,E)
    Su=quadgk(integrandseebeck,min,Ef+20kBe*Temp)[1]    
    return -1/q/Temp*(Su)
end
function seebeck_constanttau_Multiband(tau,bndst,Ef,Temp)
    seebeckTotal=0.0
    seebeckTotalDenominator=sigma_constanttau_Multiband(tau,bndst,Ef,Temp)   
    for band in bndst.bands
        seebeckTotal+=seebeck_constanttau_Nominator(tau,band,Ef,Temp)/seebeckTotalDenominator        
    end
    return seebeckTotal
end
function seebeck_Multiband(tau_electron::types.tau_electron_Base,bndst,Ef,Temp)
    seebeckTotal=0.0
    seebeckTotalDenominator=sigma_Multiband(tau_electron,bndst,Ef,Temp)
    for band in bndst.bands   
        seebeckTotal+=seebeck_Nominator(tau_electron,band,Ef,Temp)./seebeckTotalDenominator  
    end
    return seebeckTotal
end
function seebeck_MultibandFast(tau_electron::types.tau_electron_Base,bndst,Ef,Temp,sigmaM)
    seebeckTotal=0.0
    seebeckx=0.0
    bandpre=bndst.bands[1]
    seebeckTotalDenominator=sigmaM
    for (i,band) in enumerate(bndst.bands)
        if i==1
            seebeckx=seebeck_Nominator(tau_electron,band,Ef,Temp)./seebeckTotalDenominator            
        else
            seebeckx= band==bandpre ? seebeckx : seebeck_Nominator(tau_electron,band,Ef,Temp)./seebeckTotalDenominator             
        end
        seebeckTotal+=seebeckx
        bandpre=band
    end
    return seebeckTotal
end
#Seebeck
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#Thermal Conductivity

function keint(tau_electron::types.tau_electron_Base,band,Ef,Temp)
    if band.effMass>0                
        integrandp_keint(E)=sigmaD(tau_electron,band,E,Ef,Temp).*((E-Ef)*q).^2 
        nodes, weights = qnwlege(100, band.offset,band.offset+20kBe*Temp)
        a= do_quad(integrandp_keint,nodes, weights)
        #Su=a#quadgk(integrandseebeck,band.offset,band.offset+20kBe*Temp)[1]
        return a        
    elseif band.effMass<0
        min=band.offset-20kBe*Temp<0 ? 0.0 : band.offset-20kBe*Temp    
        integrandn_keint(E)=sigmaD(tau_electron,band,E,Ef,Temp).*((E-Ef)*q).^2
        nodes, weights = qnwlege(100,min,band.offset)
        a= do_quad(integrandn_keint,nodes, weights)
        #Su=a#quadgk(integrandseebeck,min,band.offset)[1]
        return a
    else
        error("Band EffMass cannot be =0")
        return 0.0
    end
    #return a#quadgk(integrand,min,Ef+20kBe*Temp)[1]
end
#
#
#
function ke(tau_electron::types.tau_electron_Base,bndst,Ef,Temp)
    Se=0.0
    Sh=0.0
    sigmae=0.0
    sigmah=0.0
    keintTotal=0.0
    for band in bndst.bands        
        if band.effMass>=0
            keintTotal+=keint(tau_electron,band,Ef,Temp)
            Se+=seebeck(tau_electron,band,bndst,Ef,Temp)
            #seebeck(tau_electron::types.tau_electron_Base,iband,bndst,Ef,Temp)
            sigmae+=sigma(tau_electron,band,Ef,Temp)
        elseif band.effMass<0
            Sh+=seebeck(tau_electron,band,bndst,Ef,Temp)
            sigmah+=sigma(tau_electron,band,Ef,Temp)
        end                    
    end
    return 1/q/q/Temp*keintTotal-Se^2*sigmae*Temp-Sh^2*sigmah*Temp    
end
#
#
#
function keint_constanttau(tau,band,Ef,Temp)
    min=Ef-10kBe*Temp<0 ? 0.0 : Ef-10kBe*Temp 
    integrand(E)=q*q*tau*(square_velocity_E(band,E))*q*getDOS_SingleBand_E(band,E)*-fermiDerivativeTemp_Ef_E(Ef,Temp,E)*((E-Ef)*q)^2        
    return quadgk(integrand,min,Ef+10kBe*Temp)[1]
end
#
#
#
function ke_constanttau(tau,bndst,Ef,Temp)
    Se=0.0
    Sh=0.0
    sigmae=0.0
    sigmah=0.0
    keintTotal=0.0
    for band in bndst.bands
        keintTotal+=keint_constanttau(tau,band,Ef,Temp)
        if band.effMass>=0
            Se+=seebeck_constanttau(tau,band,bndst,Ef,Temp)
            #Seu=seebeck_constanttau_Nominator(tau,band,Ef,Temp)
            sigmae+=sigma_constanttau(tau,band,Ef,Temp)
        elseif band.effMass<0
            Sh+=seebeck_constanttau(tau,band,bndst,Ef,Temp)
            #Shu=seebeck_constanttau_Nominator(tau,band,Ef,Temp)
            sigmah+=sigma_constanttau(tau,band,Ef,Temp)
        end                    
    end
    return 1/q/q/Temp*keintTotal-Se^2*sigmae*Temp-Sh^2*sigmah*Temp# 
    
end
#Thermal Conductivity
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#ZT

function ZT_constanttau(tau,bndst,Ef,Temp,k)
    k=ke_constanttau(tau,bndst,Ef,Temp)+k
    S=seebeck_constanttau_Multiband(tau,bndst,Ef,Temp)
    sigma=sigma_constanttau_Multiband(tau,bndst,Ef,Temp)
    return S^2*sigma/k*Temp    
end

#ZT
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#Lattice Thermal Conductivity 
#
#
function Ii1(tauPHTOT::tau_phonon_Base,thetai::Float64,T::Float64)
    tauPHTOT.variables[2]=T
    integrand(x)=get_tau_phonon(tauPHTOT,x).*(x.^4).*exp(x)./(exp(x)-1).^2
    nodes, weights = qnwlege(100,0,thetai/T)
    return a= do_quad(integrand,nodes, weights)    
end
#
#
function Ii2(tauPHTOT::tau_phonon_Base,tauPHN::tau_phonon_Base,thetai::Float64,T::Float64)
    tauPHTOT.variables[2]=T
    tauPHN.variables[2]=T
    integrand(x)=get_tau_phonon(tauPHTOT,x)./get_tau_phonon(tauPHN,x).*(x.^4).*exp(x)./(exp(x)-1).^2
    nodes, weights = qnwlege(100,0,thetai/T)
    return a= do_quad(integrand,nodes, weights)    
end
#
#
function Ii3(tauPHTOT::tau_phonon_Base,tauPHN::tau_phonon_Base,tauPHR::tau_phonon_Base,thetai::Float64,T::Float64)
    tauPHTOT.variables[2]=T
    tauPHN.variables[2]=T
    tauPHR.variables[2]=T
    integrand(x)=get_tau_phonon(tauPHTOT,x)./get_tau_phonon(tauPHN,x)./get_tau_phonon(tauPHR,x).*(x.^4).*exp(x)./(exp(x)-1).^2
    nodes, weights = qnwlege(100,0,thetai/T)
    return a= do_quad(integrand,nodes, weights)    
end
#
#
function kl(tauPHTOTL::tau_phonon_Base,tauPHNL::tau_phonon_Base,tauPHRL::tau_phonon_Base,
    tauPHTOTTx::tau_phonon_Base,tauPHNTx::tau_phonon_Base,tauPHRTx::tau_phonon_Base,
    tauPHTOTTy::tau_phonon_Base,tauPHNTy::tau_phonon_Base,tauPHRTy::tau_phonon_Base,
    T::Float64,v::Array{Float64,1})
    IL1=Ii1(tauPHTOTL,tauPHTOTL.variables[6],T)
    IL2=Ii2(tauPHTOTL,tauPHNL,tauPHTOTL.variables[6],T)
    IL3=Ii3(tauPHTOTL,tauPHNL,tauPHRL,tauPHTOTL.variables[6],T)
    ITx1=Ii1(tauPHTOTTx,tauPHTOTTx.variables[7],T)
    ITx2=Ii2(tauPHTOTTx,tauPHNTx,tauPHTOTTx.variables[7],T)
    ITx3=Ii3(tauPHTOTTx,tauPHNTx,tauPHRTx,tauPHTOTTx.variables[7],T)
    ITy1=Ii1(tauPHTOTTy,tauPHTOTTy.variables[8],T)
    ITy2=Ii2(tauPHTOTTy,tauPHNTy,tauPHTOTTy.variables[8],T)
    ITy3=Ii3(tauPHTOTTy,tauPHNTy,tauPHRTy,tauPHTOTTy.variables[8],T)
    kL=(kB^4*T^3/2/pi/pi/hbar^3)*(IL1+IL2*IL2/IL3)/v[1]
    kTx=(kB^4*T^3/2/pi/pi/hbar^3)*(ITx1+ITx2*ITx2/ITx3)/v[2]
    kTy=(kB^4*T^3/2/pi/pi/hbar^3)*(ITy1+ITy2*ITy2/ITy3)/v[3]
    kl=(kL+kTx+kTy)/3.0 
    return kl#(IL1,ITx1,IL2,ITx2,IL3,ITx3,kL/3,kTx/3,kTy/3,kl)#kl#(IL1,ITx1,IL2,ITx2,IL3,ITx3,kL/3,kTx/3,kTy/3,kl)#(kL+kTx+kTy)/3.0    
end
#Lattice Thermal Conductivity
##############################################################################################################
##############################################################################################################

##############################################################################################################
##############################################################################################################
#Temperature Depenedent Functions

function sigma_T_constanttau(tau,bndst,Ef,TempArray)
    sigma_T=Array{Float64}(length(TempArray))
    for (i,Temp) in enumerate(TempArray)
        sigma_T[i]= sigma_constanttau_Multiband(tau,bndst,Ef,Temp)
    end
    return sigma_T
end

function sigma_T(tau_electron::types.tau_electron_Base,bndst,Ef,TempArray)
    sigma_T=Array{Float64}(length(TempArray))    
    for (i,Temp) in enumerate(TempArray)
        tau_electron.variables[2]=Temp
        sigma_T[i]= sigma_Multiband(tau_electron,bndst,Ef,Temp)
    end
    return sigma_T
end

#ZT
##############################################################################################################
##############################################################################################################