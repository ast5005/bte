function sigmaArray(bndstTx,Ts,numofn,tauTOTTx,xmax)
    Efs=Array{Float64}(length(Ts),length(numofn)) 
    Efs=EfsArray(bndstTx,Ts,numofn,tauTOTTx,xmax)
    sigmas=Array{Float64}(length(Ts),length(numofn))
    for (j,Tx) in enumerate(Ts) 
        bndstTx.var[1]=Tx    
        tauTOTTx.variables[2]=Tx    
        types.updatebnstTx(bndstTx)    
        for (k,nx) in enumerate(numofn)
            tauTOTTx.variables[7]=nx*1e6
            sigmas[j,k]=sigma_Multiband(tauTOTTx,bndstTx,Efs[j,k],Tx)
            seebecks[j,k]=seebeck_MultibandFast(tauTOTTx,bndstTx,Efs[j,k],Tx,sigmas[j,k])
            #seebeck_Multiband(tauTOTTx,bndstTx,Efs[j,k],Tx)
        end
    end
    return (sigmas,seebecks)
end
function sigmaArray(bndstTx,Ts,xs,numofn,tauTOTTx,xmax)
    sigmasx=Array{Float64}(length(Ts),length(numofn),length(xs))
    seebecksx=Array{Float64}(length(Ts),length(numofn),length(xs))
    for (i,xsx) in enumerate(xs)
        bndstTx.var[2]=xsx
        tauTOTTx.variables[4]=xsx
        (sigmasx[:,:,i],seebecksx[:,:,i])=sigmaArray(bndstTx,Ts,numofn,tauTOTTx,xmax)
    end
    return (sigmasx,seebecksx)    
end
function EfsArray(bndstTx,Ts,numofn,tauTOTTx,xmax)
    Efs=Array{Float64}(length(Ts),length(numofn))    
    for (j,Tx) in enumerate(Ts) 
        bndstTx.var[1]=Tx    
        tauTOTTx.variables[2]=Tx    
        types.updatebnstTx(bndstTx)   
        for (k,nx) in enumerate(numofn)
            tauTOTTx.variables[7]=nx*1e6
            Efs[j,k]=Fermilevel_n(nx,bndstTx,Tx,xmax) # Efinput[j]+Eoff#
            #numofnsout[j,k]=NumofnMultiBand(bndstTx,Efs[j,k],Tx,xmax)    
        end
        
    end
    return Efs
end
function electronicProps(bndst,xs,Ts,numofn,tauTOTTx,xmax)
    bandpre=bndst.bands[1]
    sigmae=0.0#Array{Float64}(length(bndst.bands))
    sigmah=0.0#Array{Float64}(length(bndst.bands))
    seebecke=0.0#Array{Float64}(length(bndst.bands))
    seebeckh=0.0#Array{Float64}(length(bndst.bands))
    kee=0.0#Array{Float64}(length(bndst.bands))
    keh=0.0#Array{Float64}(length(bndst.bands))
    presigmae=0.0
    presigmah=0.0
    preseebecke=0.0
    preseebeckh=0.0
    prekee=0.0
    prekeh=0.0    
    sigmas=Array{Float64}(length(Ts),length(numofn),length(xs))
    seebecks=Array{Float64}(length(Ts),length(numofn),length(xs))
    kees=Array{Float64}(length(Ts),length(numofn),length(xs))
    kehs=Array{Float64}(length(Ts),length(numofn),length(xs))
    kes=Array{Float64}(length(Ts),length(numofn),length(xs))
    kbis=Array{Float64}(length(Ts),length(numofn),length(xs))
    Efs=Array{Float64}(length(Ts),length(numofn),length(xs))     
    for (i,Tx) in enumerate(Ts)
        bndstTx.var[1]=Tx    
        tauTOTTx.variables[2]=Tx 
        for (j,nx) in enumerate(numofn)
            tauTOTTx.variables[7]=nx*1e6
            for (k,xsx) in enumerate(xs)                
                sigmae=0.0
                sigmah=0.0
                seebecke=0.0
                seebeckh=0.0
                kee=0.0
                keh=0.0
                bndstTx.var[2]=xsx
                tauTOTTx.variables[4]=xsx
                types.updatebnstTx(bndst)
                Efs[i,j,k]=Fermilevel_n(nx,bndst,Tx,xmax)                
                for (l,band) in enumerate(bndst.bands)
                    if l==1
                        if band.effMass>0
                            presigmae=sigma(tauTOTTx,band,Efs[i,j,k],Tx)                            
                            preseebecke=seebeck_Nominator(tauTOTTx,band,Efs[i,j,k],Tx) 
                            prekee=keint(tauTOTTx,band,Efs[i,j,k],Tx) 
                            sigmae=isequal(presigmae,NaN) ? 0.0 : presigmae
                            seebecke=isequal(preseebecke ,NaN) ? 0.0 : preseebecke 
                            kee=isequal(prekee ,NaN) ? 0.0 : prekee 
                            
                        else
                            presigmah=sigma(tauTOTTx,band,Efs[i,j,k],Tx) 
                            preseebeckh=seebeck_Nominator(tauTOTTx,band,Efs[i,j,k],Tx)
                            prekeh=keint(tauTOTTx,band,Efs[i,j,k],Tx)
                            sigmah=isequal(presigmah,NaN) ? 0.0 : presigmah
                            seebeckh=isequal(preseebeckh ,NaN) ? 0.0 : preseebeckh 
                            keh=isequal(prekeh ,NaN) ? 0.0 : prekeh
                        end                        
                    else
                        if band.effMass>0
                            presigmae= bandpre==band ? presigmae : sigma(tauTOTTx,band,Efs[i,j,k],Tx)
                            preseebecke= bandpre==band ? preseebecke : seebeck_Nominator(tauTOTTx,band,Efs[i,j,k],Tx)
                            prekee= bandpre==band ? prekee : keint(tauTOTTx,band,Efs[i,j,k],Tx)
                            sigmae+=isequal(presigmae,NaN) ? 0.0 : presigmae
                            seebecke+=isequal(preseebecke ,NaN) ? 0.0 : preseebecke 
                            kee+=isequal(prekee ,NaN) ? 0.0 : prekee 
                        else
                            presigmah= bandpre==band ? presigmah : sigma(tauTOTTx,band,Efs[i,j,k],Tx)
                            preseebecke= bandpre==band ? preseebecke : seebeck_Nominator(tauTOTTx,band,Efs[i,j,k],Tx) 
                            prekeh= bandpre==band ? prekeh : keint(tauTOTTx,band,Efs[i,j,k],Tx)
                            sigmah+=isequal(presigmah,NaN) ? 0.0 : presigmah
                            seebeckh+=isequal(preseebeckh ,NaN) ? 0.0 : preseebeckh 
                            keh+=isequal(prekeh ,NaN) ? 0.0 : prekeh
                        end
                        bandpre=band
                    end
                end
                sigmas[i,j,k]=sigmae+sigmah
                seebecks[i,j,k]=(seebecke+seebeckh)/sigmas[i,j,k]            
                kbis[i,j,k]=(sigmae*sigmah)/(sigmae+sigmah)*((seebecke/sigmae-seebeckh/sigmah)^2)*Tx                
                kees[i,j,k]=kee/q/q/Tx-((seebecke/sigmae)^2)*sigmae*Tx
                kehs[i,j,k]=keh/q/q/Tx-((seebeckh/sigmah)^2)*sigmah*Tx
                kes[i,j,k]=kees[i,j,k]+kehs[i,j,k]+kbis[i,j,k]
                sigmas[i,j,k]=isequal(sigmas[i,j,k],NaN) ? 0.0 : sigmas[i,j,k]
                seebecks[i,j,k]=isequal(seebecks[i,j,k],NaN) ? 0.0 : seebecks[i,j,k]
                kbis[i,j,k]=isequal(kbis[i,j,k],NaN) ? 0.0 : kbis[i,j,k]
                kees[i,j,k]=isequal(kees[i,j,k],NaN) ? 0.0 : kees[i,j,k]
                kehs[i,j,k]=isequal(kehs[i,j,k],NaN) ? 0.0 : kehs[i,j,k]
                kes[i,j,k]=isequal(kes[i,j,k],NaN) ? 0.0 : kes[i,j,k]
            end
        end
   end
    return (Efs,sigmas,seebecks,kes,kees,kehs,kbis)
end
function electronicPropsiso(bndst,xs,Ts,numofn,tauTOTTx,xmax)
    bandpre=bndst.bands[1]
    sigmae=0.0#Array{Float64}(length(bndst.bands))
    sigmah=0.0#Array{Float64}(length(bndst.bands))
    seebecke=0.0#Array{Float64}(length(bndst.bands))
    seebeckh=0.0#Array{Float64}(length(bndst.bands))
    kee=0.0#Array{Float64}(length(bndst.bands))
    keh=0.0#Array{Float64}(length(bndst.bands))
    presigmae=0.0
    presigmah=0.0
    preseebecke=0.0
    preseebeckh=0.0
    prekee=0.0
    prekeh=0.0    
    sigmas=Array{Float64}(length(Ts),length(numofn),length(xs))
    seebecks=Array{Float64}(length(Ts),length(numofn),length(xs))
    kees=Array{Float64}(length(Ts),length(numofn),length(xs))
    kehs=Array{Float64}(length(Ts),length(numofn),length(xs))
    kes=Array{Float64}(length(Ts),length(numofn),length(xs))
    kbis=Array{Float64}(length(Ts),length(numofn),length(xs))
    Efs=Array{Float64}(length(Ts),length(numofn),length(xs))     
    for (i,Tx) in enumerate(Ts)
        bndstTx.var[1]=Tx    
        tauTOTTx.variables[2]=Tx 
        for (j,nx) in enumerate(numofn)
            k=j
            xsx=xs[j]
            tauTOTTx.variables[7]=nx*1e6
            #for (k,xsx) in enumerate(xs)                
                sigmae=0.0
                sigmah=0.0
                seebecke=0.0
                seebeckh=0.0
                kee=0.0
                keh=0.0
                bndstTx.var[2]=xsx
                tauTOTTx.variables[4]=xsx
                types.updatebnstTx(bndst)
                Efs[i,j,k]=Fermilevel_n(nx,bndst,Tx,xmax)                
                for (l,band) in enumerate(bndst.bands)
                    if l==1
                        if band.effMass>0
                            presigmae=sigma(tauTOTTx,band,Efs[i,j,k],Tx)                            
                            preseebecke=seebeck_Nominator(tauTOTTx,band,Efs[i,j,k],Tx) 
                            prekee=keint(tauTOTTx,band,Efs[i,j,k],Tx) 
                            sigmae=isequal(presigmae,NaN) ? 0.0 : presigmae
                            seebecke=isequal(preseebecke ,NaN) ? 0.0 : preseebecke 
                            kee=isequal(prekee ,NaN) ? 0.0 : prekee 
                            
                        else
                            presigmah=sigma(tauTOTTx,band,Efs[i,j,k],Tx) 
                            preseebeckh=seebeck_Nominator(tauTOTTx,band,Efs[i,j,k],Tx)
                            prekeh=keint(tauTOTTx,band,Efs[i,j,k],Tx)
                            sigmah=isequal(presigmah,NaN) ? 0.0 : presigmah
                            seebeckh=isequal(preseebeckh ,NaN) ? 0.0 : preseebeckh 
                            keh=isequal(prekeh ,NaN) ? 0.0 : prekeh
                        end                        
                    else
                        if band.effMass>0
                            presigmae= bandpre==band ? presigmae : sigma(tauTOTTx,band,Efs[i,j,k],Tx)
                            preseebecke= bandpre==band ? preseebecke : seebeck_Nominator(tauTOTTx,band,Efs[i,j,k],Tx)
                            prekee= bandpre==band ? prekee : keint(tauTOTTx,band,Efs[i,j,k],Tx)
                            sigmae+=isequal(presigmae,NaN) ? 0.0 : presigmae
                            seebecke+=isequal(preseebecke ,NaN) ? 0.0 : preseebecke 
                            kee+=isequal(prekee ,NaN) ? 0.0 : prekee 
                        else
                            presigmah= bandpre==band ? presigmah : sigma(tauTOTTx,band,Efs[i,j,k],Tx)
                            preseebecke= bandpre==band ? preseebecke : seebeck_Nominator(tauTOTTx,band,Efs[i,j,k],Tx) 
                            prekeh= bandpre==band ? prekeh : keint(tauTOTTx,band,Efs[i,j,k],Tx)
                            sigmah+=isequal(presigmah,NaN) ? 0.0 : presigmah
                            seebeckh+=isequal(preseebeckh ,NaN) ? 0.0 : preseebeckh 
                            keh+=isequal(prekeh ,NaN) ? 0.0 : prekeh
                        end
                        bandpre=band
                    end
                end
                sigmas[i,j,k]=sigmae+sigmah
                seebecks[i,j,k]=(seebecke+seebeckh)/sigmas[i,j,k]            
                kbis[i,j,k]=(sigmae*sigmah)/(sigmae+sigmah)*((seebecke/sigmae-seebeckh/sigmah)^2)*Tx                
                kees[i,j,k]=kee/q/q/Tx-((seebecke/sigmae)^2)*sigmae*Tx
                kehs[i,j,k]=keh/q/q/Tx-((seebeckh/sigmah)^2)*sigmah*Tx
                kes[i,j,k]=kees[i,j,k]+kehs[i,j,k]+kbis[i,j,k]
                sigmas[i,j,k]=isequal(sigmas[i,j,k],NaN) ? 0.0 : sigmas[i,j,k]
                seebecks[i,j,k]=isequal(seebecks[i,j,k],NaN) ? 0.0 : seebecks[i,j,k]
                kbis[i,j,k]=isequal(kbis[i,j,k],NaN) ? 0.0 : kbis[i,j,k]
                kees[i,j,k]=isequal(kees[i,j,k],NaN) ? 0.0 : kees[i,j,k]
                kehs[i,j,k]=isequal(kehs[i,j,k],NaN) ? 0.0 : kehs[i,j,k]
                kes[i,j,k]=isequal(kes[i,j,k],NaN) ? 0.0 : kes[i,j,k]
            #end
        end
   end
    return (Efs,sigmas,seebecks,kes,kees,kehs,kbis)
end
function thermalPropsiso(tauPHL::Array{tau_phonon_B,1},tauPHTx::Array{tau_phonon_B,1},tauPHTy::Array{tau_phonon_B,1},vsos::Array{Float64,1})
    
    klattice=Array{Float64,3}(length(Ts),length(numofn),length(xs))
    for (i,iT) in enumerate(Ts)
    #for (k,kx) in enumerate(xs) 
        for (j,jn) in enumerate(numofn)
            k=j
            kx=xs[j]
            tauPHTOTL.variables[2]=iT 
            tauPHTOTTx.variables[2]=iT
            tauPHTOTTy.variables[2]=iT
            tauPHTOTL.variables[5]=kx 
            tauPHTOTTx.variables[5]=kx
            tauPHTOTTy.variables[5]=kx
            tauPHTOTL.variables[2]=iT 
            tauPHTOTTx.variables[2]=iT
            tauPHTOTTy.variables[2]=iT
            tauPHTOTL.variables[26]=j 
            tauPHTOTTx.variables[26]=j
            tauPHTOTTy.variables[26]=j
            tauPHTOTL.variables[27]=k 
            tauPHTOTTx.variables[27]=k
            tauPHTOTTy.variables[27]=k
            #klattice[i,k,j]=kl(tauPHTOTL,tauPHNL,tauPHRL,tauPHTOTTx,tauPHNTx,tauPHRTx,tauPHTOTTy,tauPHNTy,tauPHRTy,iT,
        #[vSiGeL_PH(xsp),vSiGeT_PH(xsp),vSiGeT_PH(xsp)])
            klattice[i,j,k]=kl(tauPHL[1],tauPHL[2],tauPHL[3],tauPHTx[1],tauPHTx[2],tauPHTx[3],tauPHTy[1],tauPHTy[2],
            tauPHTy[3],iT,vsos)
        end
    #end
    end
    return klattice#(Efs,sigmas,seebecks,kes,kees,kehs,kbis,klattice)#
end

function electronicthermalPropsiso(bndstTx,xs,Ts,numofn,tauTOTTx,xmax,
    tauPHL::Array{tau_phonon_B,1},tauPHTx::Array{tau_phonon_B,1},tauPHTy::Array{tau_phonon_B,1},vsos::Array{Float64,1})
    (Efs,sigmas,seebecks,kes,kees,kehs,kbis)=electronicPropsiso(bndst,xs,Ts,numofn,tauTOTTx,xmax)
klattice=Array{Float64,3}(length(Ts),length(numofn),length(xs))
for (i,iT) in enumerate(Ts)
    #for (k,kx) in enumerate(xs) 
        for (j,jn) in enumerate(numofn)
            k=j
            kx=xs[j]
            tauPHTOTL.variables[2]=iT 
            tauPHTOTTx.variables[2]=iT
            tauPHTOTTy.variables[2]=iT
            tauPHTOTL.variables[5]=kx 
            tauPHTOTTx.variables[5]=kx
            tauPHTOTTy.variables[5]=kx
            tauPHTOTL.variables[2]=iT 
            tauPHTOTTx.variables[2]=iT
            tauPHTOTTy.variables[2]=iT
            tauPHTOTL.variables[26]=j 
            tauPHTOTTx.variables[26]=j
            tauPHTOTTy.variables[26]=j
            tauPHTOTL.variables[27]=k 
            tauPHTOTTx.variables[27]=k
            tauPHTOTTy.variables[27]=k
            #klattice[i,k,j]=kl(tauPHTOTL,tauPHNL,tauPHRL,tauPHTOTTx,tauPHNTx,tauPHRTx,tauPHTOTTy,tauPHNTy,tauPHRTy,iT,
        #[vSiGeL_PH(xsp),vSiGeT_PH(xsp),vSiGeT_PH(xsp)])
            klattice[i,j,k]=kl(tauPHL[1],tauPHL[2],tauPHL[3],tauPHTx[1],tauPHTx[2],tauPHTx[3],tauPHTy[1],tauPHTy[2],
            tauPHTy[3],iT,vsos)
        end
    #end
end
    return (Efs,sigmas,seebecks,kes,kees,kehs,kbis,klattice)#klattice
end
