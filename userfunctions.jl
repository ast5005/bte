#userfunctions.jl

using constants
using Cubature
using QuantEcon

function gen(x)
    return 1
end

global g_T=0.0
global g_Ef=0.0
global g_band=parBandTx(0.0,0.0,0.0,[func1(x)=1],[func2(x)=1],[0.0,0.0])
global g_intvalue=0.0
#global g_Ef=0.0,g_band=parBandTx(0.0,0.0,0.0,[gen],[gen],[0.0,0.0]),g_intvalue=0.0
#global g_Ef=0.0,g_band=parBandTx(0.0,0.0,0.0,[function ft1(x) return 1 end],[function ft2(x) return 1 end],[0.0,0.0]),g_intvalue=0.0

function Egx3(T,x)
    return (0.78-(4.0e-4)*T)*(1-x)+(0.38-(2.8e-4)*T)*x
end

function Ec0(x)
    return 0.4*(1-x)-0.165*x 
end

function X1effMass(T)
    return 0.49+2.0e-4*T 
end

function epsilon0_x(x)
    return (20*(1-x)+23.75*x)*8.85e-12 
end

function epsilonhf_x(x)
    return (13.3*(1-x)+17*x)*8.85e-12
end
function opPhE(x)
    return (40.5*(1-x)+28.8*x)/1000 
end
function LatticeConst(x)
    return 6.338*(1-x)+3.59e3*x   
end
function acPhDa(x)
    return 9.0-1.5*x*x 
end
function Cl(x)
    return 4.15e10*(1-x)+3.22e10*x 
end
######################################################
##SiGe Properties
function epsilon_SiGe()
    return 27.4; 
end
function EgSi(T::Float64)
    return (1.17-(4.73e-4*(T*T))/(T+636))
end
function EgGe(T::Float64)
    return (0.742-(4.8e-4*(T*T))/(T+235))
end
function EgSiGe(T::Float64,x::Float64)
    return EgSi(T)*(1-x)+EgGe(T)*x
    #(0.8941+0.042*(1-x)+0.1961*(1-x)^2-(0.00037*x+0.00023*(1-x)));
end
function ClSiGe(x::Float64,T::Float64)
    return vSiGe(x).^2.*roSiGe(x)
end
function DefP(md::Float64)
    return  15#7.0
end
function VSiGe(x::Float64)
    return ((2.7155e-10)^3*(1-x)+(2.8288e-10)^3*x)
end
function MSiGe(x::Float64)    
    return (28.086e-3*(1-x)+72.59e-3*x)/Nav
end
function vSiGe(x::Float64)
    #return kB*((VSiGe(x)./6/pi/pi)^(1.0/3.0))*586/hbar
    return kB*((VSiGe(x)./6/pi/pi)^(1.0/3.0))*thetaAlloySiGe(x)/hbar
end
function thetaAlloySiGe(x::Float64)
    Nav=6.0221415e23
    G=1.033*(1-x)+1.017*x    
    a=VSiGe(x)^(1/3)    
    return 1.48e-11*a^(-3/2)*(MSiGe(x)*Nav*1000)^(-1/2)*G
end
function roSiGe(x::Float64)
    return MSiGe(x)/VSiGe(x)
end
function C11SiGe(x::Float64)
    return 165.8-37.3x
end
function C44SiGe(x::Float64)
    return 79.6-12.8x
end
function vSiGeL(x::Float64)
    roSiGeSI=(2.329+3.493*x-0.499*x*x)
    return 700*(C11SiGe(x)/roSiGeSI)^0.5    
end
function vSiGeT(x::Float64)
    roSiGeSI=(2.329+3.493*x-0.499*x*x)
    return 700*(C44SiGe(x)/roSiGeSI)^0.5    
end
function vSiGeL_PH(x::Float64)
    return 8430*(1-x)+4920*x
end
function vSiGeT_PH(x::Float64)
    return 5840*(1-x)+3540*x
end
function thetaCL(a::Float64,x::Float64,thetaD::Float64,vSi::Float64,VSi::Float64)
    return a*thetaD*vSi/vSiGeL_PH(x)*(VSi/VSiGe(x))^1/3
end
function thetaCT(a::Float64,x::Float64,thetaD::Float64,vSi::Float64,VSi::Float64)
    return a*thetaD*vSi/vSiGeT_PH(x)*(VSi/VSiGe(x))^1/3
end
function Gamma(MI::Float64,MII::Float64,aI::Float64,aII::Float64,xs::Float64,eps::Float64)
    delM=MI-MII
    M=MI*(1-xs)+MII*xs
    a=aI*(1-xs)+aII*xs
    dela=aI-aII
   return xs*(1-xs)*((delM/M)^2+eps*(dela/a)^2) 
end
##SiGe Properties
######################################################
function tauAC_func(Cl,Da,T,mdx,E::Float64,band)    
    md=band.effMass
    t=1#1e1
    if md>0   
        #println("tauAC_func md>0")
        if E-band.offset < 0 
            #println("tauAC_func E-band.offset<0")
            return 0    
        else        
            #println("tauAC_func E-band.offset>0")      
            
            return hbar^4*Cl*pi/(sqrt(2*(E-band.offset)*q).*abs(md)^1.5.*(Da*q)^2.*(kB*T))
        end
    else
        #println("tauAC_func md<0")
        if -E+band.offset < 0 
            #println("tauAC_func -E+band.offset<0")
            #println("E= $E")
            #println("Offset=$(band.offset)")            
            return 0    
        else        
            #println("tauAC_func -E+band.offset>0")
            return hbar^4*Cl*pi/(sqrt(2*(-E+band.offset)*q)*abs(md)^1.5*(Da*q)^2*(kB*T))
        end
    end
end
function tauAC_func(Cl::Float64,Da::Float64,Dv::Float64,T::Float64,mdx::Float64,E::Float64,band::parBandTx)      
    a=band.alpha
    md=band.effMass
    t=1
    Da=Da*q
    Dv=Dv*q
    if md>0  
        if E-band.offset < 0 
            return 0    
        else      
            Et=(E-band.offset)*q     
            aq=a/q
            Eta=Et*aq
            De=(1-Eta/(1+2*Eta)*(1-Dv/Da))^2-8/3*(Eta*(1+Eta)*Dv)/((1+2*Eta)^2*Da)            
            return 1./(pi*kB*T*Da^2/Cl/hbar*getDOS_SingleBand_E(band,E)*De) 
            #hbar^4*Cl*pi/(sqrt(2*((E-band.offset)*q+aq*Et^2))*(1+2*Et*aq).*abs(md)^1.5*(Da*q)^2.*(kB*T))/De
        end
    else
        if -E+band.offset < 0
            return 0    
        else
            Et=(-E+band.offset)*q
            aq=a/q
            Eta=Et*a
            De=(1-Eta./(1+2*Eta)*(1-Dv/Da))^2-8/3*(Eta*(1+Eta)*Dv)/((1+2*Eta)^2*Da)
            return 1./(pi*kB*T*Da^2/Cl/hbar*getDOS_SingleBand_E(band,E)*De) 
            #hbar^4*Cl*pi/(sqrt(2*((-E+band.offset)*q+aq*Et^2))*(1+2*Et*aq)*abs(md)^1.5*(Da*q)^2*(kB*T))/De
        end
    end
end
function tauAC_func(Cl,T,md,E::Float64,band)  
    #println("Reached tauAC_func")
    t=1#1e1
    #foreach(println,band.var)
    Da=band.var[3]
    if md>0   
        #println("tauAC_func md>0")
        if E-band.offset < 0 
            #println("tauAC_func E-band.offset<0")
            return 0    
        else        
            #println("tauAC_func E-band.offset>0")           
            return hbar^4*Cl*pi/(sqrt(2*(E-band.offset)*q)*abs(md)^1.5*(Da*q)^2*(kB*T))
        end
    else
        #println("tauAC_func md<0")
        if -E+band.offset < 0 
            #println("tauAC_func -E+band.offset<0")
            #println("E= $E")
            #println("Offset=$(band.offset)")            
            return 0    
        else        
            #println("tauAC_func -E+band.offset>0")
            return hbar^4*Cl*pi/(sqrt(2*(-E+band.offset)*q)*abs(md)^1.5*(Da*q)^2*(kB*T))
        end
    end
end
function tauAC_func(Cl,Da,T,md,E::Array{Float64},band)
    tau=Array{Float64}(length(E))
    for (i,Ex) in enumerate(E)
        temp= tauAC_func(Cl,Da,T,md,Ex,band) 
        tau[i]=  temp>0 ? temp :0.0 
    end
    return tau
end
function tauAC_func(Cl::Float64,Da::Float64,Dv::Float64,T::Float64,md::Float64,E::Array{Float64},band::parBandTx)
    tau=Array{Float64}(length(E))
    for (i,Ex) in enumerate(E)
        temp= tauAC_func(Cl,Da,Dv,T,md,Ex,band) 
        tau[i]=  temp>0.0 ? temp :0.0 
    end
    return tau
end
function tauAC_func(Cl,T,md,E::Array{Float64},band)
    tau=Array{Float64}(length(E))
    for (i,Ex) in enumerate(E)
        temp= tauAC_func(Cl,T,md,Ex,band) 
        tau[i]=  temp>0 ? temp :0.0 
    end
    return tau
end
function tauPOP_func(epsilon0,epsilonhf,band,Ef,T,md,E,integrali,NII,opPhE,bndst)
    t=1#0.05#1.0e-77
    N0=1/(exp(q*opPhE/kB/(T))-1.0)    
    min=Ef-12kBe*T<0 ? 0.0 : Ef-12kBe*T
    integrand(x)=x<band.offset ? 0.0 : q*getDOS_MultiBand_E_Total(bndst,x)*-fermiDerivativeTemp_Ef_E(Ef,T,x) 
    integral=quadgk(integrand,min,(Ef+20kBe*T))[1]
    rhfsq=epsilonhf/(q*q*integral)#epsilon0*kB*T/q/q/abs(NII)#2*epsilon0*(Ef-band.offset)*q/(3*q*q*abs(NII))
    if md>0
      if E-band.offset < 0 
            return 0    
        else 
            Ex=E-band.offset
            deltahf=hbar*hbar/(8*md*rhfsq*(E-band.offset)*q)    
            brac=Ex>opPhE ? N0*sqrt(1+opPhE/Ex)+(N0+1)*(sqrt(1-opPhE/Ex))-((N0)*opPhE/Ex)*asinh((Ex/opPhE)^0.5)+((N0+1)*opPhE/Ex)*asinh((Ex/opPhE-1)^0.5) : 1
            return 10/(1-deltahf*log(1+1/deltahf))*hbar*hbar*sqrt((E-band.offset)*q)/(sqrt(2*md)*q*q*opPhE*q)/((1.0/epsilonhf)-(1.0/epsilon0))#t*hbar*hbar*sqrt((E-band.offset)*q)/(sqrt(2*md)*q*q*kB*T*((1.0/epsilonhf)-(1.0/epsilon0)))#/(1-deltahf*log(1+1/deltahf))
        end
    else
        if -E+band.offset < 0 
            return 0 
        else        
            deltahf=hbar*hbar/(8*md*rhfsq*(-E+band.offset)*q)          
            return t*hbar^2*sqrt((-E+band.offset)*q/2/md)/(q*q*kB*T*(1/epsilonhf-1/epsilon0))/(1-deltahf*log(1+1/deltahf))
        end
    end
end

function tauII_func(epsilon0,epsilonhf,band,Ef,T,NII,md,E,integral)
    t=1#0.1e-18#7e-2
    #println("buraya geldi")
    if md>0
        if E-band.offset < 0 
            #println("maalesef burada")
            return 0 
            
        else              
            r0sq=1/((q*q/epsilon0)*integral) #1/(q*q*abs(NII)/kB/T/epsilon0) #          
            delta0=(hbar^2)/(8*md*r0sq*(E-band.offset)*q)            
            return t*16*sqrt(2*md)*pi*epsilon0^2*((E-band.offset)*q)^1.5/abs(NII)/q^4/(log(1+1/delta0)-delta0/(1+delta0))#delta0/(1+delta0)
        end
    else
        if -E+band.offset < 0 
            #println("-E+band.offset < 0")
            return 0    
        else             
            r0sq=1/((q*q/epsilon0)*integral) #1/(q*q*abs(NII)/kB/T/epsilon0)# q*q*abs(NII)/kB/T/epsilon0#   1/((q*q/epsilon0)*integral)#         
            delta0=(hbar^2)/(8*md*r0sq*(-E+band.offset)*q)           
            #println("Here")
            return t*16*sqrt(2*md)*pi*epsilon0^2*((-E+band.offset)*q)^1.5/abs(NII)/q^4/(log(1+1/delta0)-delta0/(1+delta0))
        end
    end    
end

function tauPOPIIint(Ef,T,x,band)
    min=Ef-12kBe*T<0 ? 0.0 : Ef-12kBe*T
    integrand(x)=q*getDOS_SingleBand_E(band,x)*-fermiDerivativeTemp_Ef_E(Ef,T,x)  
    return quadgk(integrand,min,(Ef+20kBe*T))[1]#0.75*(T-300)*1e45/600+0.5e45 # hquadrature(integrand, min, Ef+12kBe*T, reltol=1e-4, abstol=1e-4, maxevals=0)[1]#quadgk(integrand,min,(Ef+8kBe*T))[1]#   
end

function tauPOP2_func(epsilon0,epsilonhf,band,Ef,T,md,E::Float64,opPhE)    
    t=1
    #opPhE=opPhE
    Ex=0
    brac=0
    epsilon=8.85e-12
    N0=1/(exp(q*opPhE/kB/(T))-1.0)      
    mul=(4*pi*epsilon0*hbar)/(q*q*q*opPhE/hbar*(epsilon0/epsilonhf-1))#(4*pi*epsilon0*hbar)/(q*q*q*opPhE/hbar*((epsilon0/epsilonhf)-1))
    if md>0
        if E-band.offset < 0 
            return 0    
        else     
            mul=mul*sqrt(2*(E-band.offset)*q/md)
            Ex=E-band.offset 
        end    
    elseif md<0
        if -E+band.offset < 0 
            return 0    
        else     
            mul=mul*sqrt(-2*(-E+band.offset)*q/md)
            Ex=-E+band.offset
        end
    else
        mul=1
    end
    brac=Ex>opPhE ? N0*sqrt(1+opPhE/Ex)+(N0+1)*(sqrt(1-opPhE/Ex))-((N0)*opPhE/Ex)*asinh((Ex/opPhE)^0.5)+((N0+1)*opPhE/Ex)*asinh((Ex/opPhE-1)^0.5) : 1    
    return mul/brac
end

function tauPOP2_func(epsilon0,epsilonhf,band,Ef,T,md,E::Array{Float64},opPhE)
    tau=Array{Float64}(length(E))
    for (i,Ex) in enumerate(E)
        temp= tauPOP2_func(epsilon0,epsilonhf,band,Ef,T,md,Ex,opPhE) 
        tau[i]=  temp>0 ? temp :0.0 
    end
    return tau
end

function tauII2_func(epsilon0,epsilonhf,band,Ef,T,NII,md,E::Float64,bndst)
    global g_T,g_Ef,g_band,g_intvalue    
    t=1#0.1e-18#7e-2
    min=Ef-20kBe*T<0 ? 0.0 : Ef-20kBe*T    
   #epsilon0*kB*T/q/q/abs(NII)#2*epsilon0*(Ef-band.offset)*q/(3*q*q*abs(NII))
    if md>0
        if E-band.offset < 0 
            return 0    
        else  
            integrandp_tauII2_func(x)=q*getDOS_MultiBand_E_Total(bndst,x).*-fermiDerivativeTemp_Ef_E(Ef,T,x)    
            nodes, weights = qnwlege(80, min, (Ef+30kBe*T))
            #integral1=quadgk(integrand,min,Ef+20kBe*T)[1]            
            if band==g_band && Ef==g_Ef && g_T==T                
                integral2=g_intvalue               
                #println("hereC1")
            else                
                g_intvalue=do_quad(integrandp_tauII2_func,nodes, weights)
                integral2=g_intvalue
                g_band=band
                g_Ef=Ef
                g_T=T
                #println("hereC2")
            end
            r0sq=epsilon0/(q*q*integral2)
            #println("integral1=$integral1")
            #println("integral2=$integral2")
            gamma=(8*md*r0sq*(E-band.offset)*q)/(hbar^2)
            return t*16*sqrt(2*md)*pi*epsilon0^2*((E-band.offset)*q)^1.5/abs(NII)/q^4/(log(1+gamma)-gamma/(1+gamma))#delta0/(1+delta0)
            #return t*16*sqrt(2*md)*pi*epsilon0^2*((E-band.offset)*q)^1.5/abs(NII)/q^4/(log(1+gamma)-1/(1+gamma))#delta0/(1+delta0)
            #return t*16*sqrt(2*md)*pi*epsilon0^2*((E-band.offset)*q)^1.5/abs(NII)/q^4/1 #delta0/(1+delta0)
            
        end
    else
        if -E+band.offset < 0 
            return 0    
        else   
            integrandn_tauII2_func(x)=q*getDOS_MultiBand_E_Total(bndst,x).*-fermiDerivativeTemp_Ef_E(Ef,T,x)    
    nodes, weights = qnwlege(60, min, (Ef+20kBe*T))
            #integral1=quadgk(integrand,min,Ef+20kBe*T)[1]
            global g_T,g_Ef,g_bndst,g_intvalue
            if  band==g_band && Ef==g_Ef && g_T==T                
                integral2=g_intvalue
                #println("hereV1")
            else                
                g_intvalue=do_quad(integrandn_tauII2_func,nodes, weights)
                integral2=g_intvalue
                g_band=band 
                g_Ef=Ef
                g_T=T                
                #println("hereV2")
            end
            #println("integral1=$integral1")
            #println("integral2=$integral2")
            r0sq=epsilon0/(q*q*integral2)
            gamma=(8*abs(md)*r0sq*(-E+band.offset)*q)/(hbar^2)
            return t*16*sqrt(2*abs(md))*pi*epsilon0^2*((-E+band.offset)*q)^1.5/abs(NII)/q^4/(log(1+gamma)-gamma/(1+gamma))
            #return t*16*sqrt(2*abs(md))*pi*epsilon0^2*((-E+band.offset)*q)^1.5/abs(NII)/q^4/1
        end
    end    
    
end

function tauII2_func(epsilon0,epsilonhf,band,Ef,T,NII,md,E::Array{Float64},bndst)
    tau=Array{Float64}(length(E))
    temp=0.0
    for (i,Ex) in enumerate(E)
        temp=tauII2_func(epsilon0,epsilonhf,band,Ef,T,NII,md,Ex,bndst)
        tau[i]=  temp>0 ? temp :0.0       
    end
    return tau
end
function tauNI_func(epsilon0::Float64,md::Float64,NNI::Float64)
    return (md*q)^2/(20*hbar^3*(4*pi*epsilon0)*NNI)
end
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#Phonon Relaxation Times
function tauPH_UV(gamma::Float64,x::Array{Float64},T::Float64,M::Float64,theta::Float64,a::Float64,beta::Float64)
    bb=20*pi/3*hbar*(6*pi*pi/4)^(1/3)        
    #println("bb ",bb)
    return 1./(bb*(1+beta*5/9)/(1+beta)*gamma^2/M/a^2*(T/theta)^3.*x.^2)    
end
function tauPH_PDV(GM::Float64,T::Float64,x::Array{Float64},v::Float64,delta::Float64)    
    return 1./(GM/4/pi*(delta/v)^3*(kB*T/hbar)^4.*x.^4)
end
function tauPH_eV(bandst::BandStrucTx,Eep::Float64,md::Float64,x::Array{Float64},Ef::Float64,v::Float64,T::Float64,d::Float64)
    tau=0
    for band in bandst.bands
        if band.effMass>0     
            md=band.effMass
            #println("md ",md)
            ts=md*v*v/(2*kB*T)
            Efst=Ef#*q/kB/T
            ex1=(1+exp(-ts+Efst-x.*x/16/ts+x/2))
            ex2=(1+exp(-ts+Efst-x.*x/16/ts-x/2)) 
            kk=((Eep^2*md^3*v)/(4*pi*hbar^4*d*ts).*log(ex1./ex2))
            tau+=kk
            #println("kk ",1./kk)
        else
            tau+=0
        end
    end
    return 1./tau#1./((Eep^2*md^3*v)/(4*pi*hbar^4*d*ts).*log(ex1./ex2))
end
function tauPH_NL(gamma::Float64,V::Float64,x::Array{Float64},T::Float64,M::Float64,v::Float64)
    
    return M*v^5*hbar^4./(kB^5*gamma^2*V*x.*x*T^5)    
end
function tauPH_NT(gamma::Float64,V::Float64,x::Array{Float64},T::Float64,M::Float64,v::Float64)
    
    return M*v^5*hbar^4./(kB^5*gamma^2*V*x*T^5)    
end
function tauPH_U(gamma::Float64,x::Array{Float64},T::Float64,M::Float64,v::Float64,theta::Float64)
  
    return (M*v^2*theta*hbar)./((gamma*kB*x.*T).^2*T.*exp(-theta/3./T))
end
function tauPH_U_SA(gamma::Float64,intx::Array{Float64},T::Float64,M::Float64,theta::Float64,omegaD::Float64,beta::Float64,delta::Float64)#intx is different from x in other phonon calculations intx=omega/omegaD
    Mcgs=M*1000
    deltacgs=delta*100
    return 3.264e-2*((1+5/9*beta)/(1+beta))*gamma^2.*intx.^2./Mcgs/deltacgs^2/(theta/T)    
end
function tauPH_EP_SA(Eep::Float64,md::Float64,x::Array{Float64},Ef::Float64,
    T::Float64,M::Float64,theta::Float64,delta::Float64)
    eta=Ef#Ef*q/kB/T
    deltacgs=delta*100
    Mcgs=M
    A=6.76e26(md/me)^2*deltacgs^2/Mcgs
    y=3.72e9(md/me)*deltacgs^2*theta
    D=1.68e-11/(md/me)/deltacgs^2/theta
    alphat=thetaD/T
    #println("alphat ",alphat)
    lambda=1
    ex1=1+exp(-alphat*y+eta-D*alphat*x.*x+alphat.*x/2)
    ex2=1+exp(-alphat*y+eta-D*alphat*x.*x-alphat.*x/2)    
    ext=ex1./ex2
    logt=log(ext)
    return 1./(lambda*(A*Eep^2/alphat).*logt) 
end
function tauPH_ALL(MI::Float64,MII::Float64,y::Float64,V::Float64,x::Array{Float64},T::Float64,M::Float64,v::Float64)
     Gamma=(1-y)*((MI-M)/M)^2+y*((MII-M)/M)^2
    return 1./(Gamma*V/v^3./4./pi)./(kB*T*x/hbar).^4 #4*pi*v^3*hbar^4./(Gamma*V*(kB*T*x).^4)
end
function tauPH_B(v::Float64,alpha::Float64,d::Float64)
    return d*(1+alpha)/(v*(1-alpha)) 
end
function tauPH_C(Nc::Float64,Dc::Float64,x::Array{Float64},T::Float64,v::Float64)
    return (2*(kB*Dc*x*T).^4+32*hbar^4*v^4)./(pi*v*Nc*Dc^2*(Dc*kB*x*T).^4) 
end
function tauPH_e(Eep::Float64,md::Float64,ro::Float64,Ef::Float64,x::Array{Float64},T::Float64,v::Float64)
    Ef=Ef*q/kB/T
    beta=md*v^2/2/kB/T
    logp=log.((1+exp.(beta-Ef+x.*x/(16*beta)+x/2))./(1+exp.(beta-Ef+x.*x/(16*beta)-x/2)))
    return 4*pi*hbar^4*ro*beta./(3*Eep^2*md^3*v)./(x-logp) 
end
