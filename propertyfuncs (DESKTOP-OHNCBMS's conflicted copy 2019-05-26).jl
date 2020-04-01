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
function fitsigmas(bndst,xs,TsA,numofn,numofnn,tauTOTTx,xmax,sigmasA)
    bandpre=bndst.bands[1]
    sigmae=0.0#Array{Float64}(length(bndst.bands))
    sigmah=0.0#Array{Float64}(length(bndst.bands))
    presigmae=0.0
    presigmah=0.0
    numofns=Array{Float64}(length(TsA),length(numofn),length(xs))   
    sigmas=Array{Float64}(length(TsA),length(numofn),length(xs))   
    sigmaes=Array{Float64}(length(Ts),length(numofn),length(xs))
    sigmahs=Array{Float64}(length(Ts),length(numofn),length(xs))    
    Efs=Array{Float64}(length(TsA),length(numofn),length(xs))
    delErr=1.05    
    for (i,Tx) in enumerate(TsA)   
        bndstTx.var[1]=Tx    
        tauTOTTx.variables[2]=Tx 
        for (j,nxx) in enumerate(numofn)   
            Err=0.0
            maxiter=10
            iter=0
            switch=true
            nx=nxx
            k=j
            while(switch && iter<maxiter)                
                iter=iter+1
                nx=Err==0 ? nxx : (nxx*(Err)) 
                numofns[i,j,k]=nx
                #println("Error=$Err, n=$nx")                
                xsx=xs[j]
                tauTOTTx.variables[7]=nx*1e6
                tauTOTTx.variables[9]=numofnn[j]*1e6
                #for (k,xsx) in enumerate(xs)                
                sigmae=0.0
                sigmah=0.0
                bndstTx.var[2]=xsx
                #println("x=$xsx")
                tauTOTTx.variables[4]=xsx
                types.updatebnstTx(bndst)
                #println("i=$i, j=$j, k=$k")
                Efs[i,j,k]=Fermilevel_n(nx,bndst,Tx,xmax)   
                for (l,band) in enumerate(bndst.bands)
                     if l==1 
                        if band.effMass>0
                            presigmae=sigma(tauTOTTx,band,Efs[i,j,k],Tx)                                                        
                            sigmae=isequal(presigmae,NaN) ? 0.0 : presigmae                           
                        else
                            presigmah=sigma(tauTOTTx,band,Efs[i,j,k],Tx)                            
                            sigmah=isequal(presigmah,NaN) ? 0.0 : presigmah                            
                        end                        
                    else
                        if band.effMass>0
                            presigmae= bandpre==band ? presigmae : sigma(tauTOTTx,band,Efs[i,j,k],Tx)                            
                            sigmae+=isequal(presigmae,NaN) ? 0.0 : presigmae                            
                        else
                            presigmah= bandpre==band ? presigmah : sigma(tauTOTTx,band,Efs[i,j,k],Tx)                           
                            sigmah+=isequal(presigmah,NaN) ? 0.0 : presigmah
                        end
                    end
                        bandpre=band
                end              
                sigmas[i,j,k]=sigmae+sigmah                
                sigmas[i,j,k]=isequal(sigmas[i,j,k],NaN) ? 0.0 : sigmas[i,j,k] 
                Err=sigmasA[i,j,k]/sigmas[i,j,k]   
                #println("iter=$iter")    
                #println("sigmasA ",sigmasA[i,j,k], "  sigmas ",sigmas[i,j,k]  )    
                switch=(Err>delErr || Err<(1/delErr))       
                #println(" at the end Error=$Err, n=$nx , switch=$switch")
            end            
        end        
    end    
    return (sigmas,numofns)
end

function electronicPropsisoVarAll(bndst,xs,Ts,TsA,sigmasA,numofn0,numofnn,tauTOTTx,xmax)
    sigmas=Array{Float64}(length(Ts),length(numofn0),length(xs))
    sigmaes=Array{Float64}(length(Ts),length(numofn0),length(xs))
    sigmahs=Array{Float64}(length(Ts),length(numofn0),length(xs))
    seebecks=Array{Float64}(length(Ts),length(numofn0),length(xs))
    seebeckes=Array{Float64}(length(Ts),length(numofn0),length(xs))
    seebeckhs=Array{Float64}(length(Ts),length(numofn0),length(xs))
    kees=Array{Float64}(length(Ts),length(numofn0),length(xs))
    kehs=Array{Float64}(length(Ts),length(numofn0),length(xs))
    kes=Array{Float64}(length(Ts),length(numofn0),length(xs))
    kbis=Array{Float64}(length(Ts),length(numofn0),length(xs))
    Efs=Array{Float64}(length(Ts),length(numofn0),length(xs))  
    numofn=Array{Float64}(length(Ts),length(xs))  
    for (i,xsx) in enumerate(xs)
        (sigma1,numofn1)=fitsigmas(bndstTx,xsx,TsA[:,i],numofn0[i],numofnn,tauTOTTx,xmax,sigmasA[:,i])
        (numofn[:,i],Efs[:,i,i],sigmas[:,i,i],seebecks[:,i,i],kes[:,i,i],kees[:,i,i],kehs[:,i,i],kbis[:,i,i],sigmaes[:,i,i],sigmahs[:,i,i],seebeckes[:,i,i],seebeckhs[:,i,i])=electronicPropsisoVarn(bndst,xs[i],Ts,TsA[:,i],numofn1,numofnn[i],tauTOTTx,xmax)
    end
    return (numofn,Efs,sigmas,seebecks,kes,kees,kehs,kbis,sigmaes,sigmahs,seebeckes,seebeckhs)
end

function electronicPropsisoVarn(bndst,xs,Ts,TsA,numofnA,numofnn,tauTOTTx,xmax)
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
    sigmas=Array{Float64}(length(Ts))
    seebecks=Array{Float64}(length(Ts))
    sigmaes=Array{Float64}(length(Ts))
    sigmahs=Array{Float64}(length(Ts))
    seebeckes=Array{Float64}(length(Ts))
    seebeckhs=Array{Float64}(length(Ts))            
    kees=Array{Float64}(length(Ts))
    kehs=Array{Float64}(length(Ts))
    kes=Array{Float64}(length(Ts))
    kbis=Array{Float64}(length(Ts))
    Efs=Array{Float64}(length(Ts)) 
    numofn=Array{Float64}(length(Ts))  
    index=1
    for (i,Tx) in enumerate(Ts)
        #println("#########################")
        #println("Tx= ",Tx)
        #println("index= ",index)        
        previous=index
        next=1+index        
        for (iTsA,TsAx) in enumerate(TsA)
            #println("++++++++++++++++++++")            
            #println("TsAx ",TsAx)
            if iTsA==1 
                #println("iTsA==1 start")
                if Tx<=TsAx 
                    numofn[index]=numofnA[1]
                    #println("Tx<=TsAx and numofn",numofn[index])
                    index=index+1
                    break             
                elseif Tx>TsAx && Tx<TsA[2]
                   numofn[index]=(numofnA[2]-numofnA[1])/(TsA[2]-TsAx)*(Tx-TsAx)+numofnA[1]                   
                   #println("Tx>TsAx and Tx<TsA[2] and numofn",numofn[index])
                   index=index+1
                   break                
                else 
                    #println("Tx>TsAx and Tx>TsA[2] also")
                end
                #println("iTsA==1 end")
            elseif iTsA==length(TsA)
                #println("iTsA==length(TsA) started")
                if Tx>=TsAx
                   numofn[index]=numofnA[length(TsA)]                   
                   #println("Tx>=TsAx and numofn",numofn[index])
                   index=index+1
                   break               
               else
                   last=length(TsA)
                   numofn[index]=(numofnA[last]-numofnA[last-1])/(TsA[last]-TsA[last-1])*(Tx-TsA[last-1])+numofnA[last-1]                   
                   #println("Tx<TsAx and numofn",numofn[index])
                   index=index+1
                   break 
               end 
                #println("iTsA==length(TsA) ended")
            else
                #println("iTsA between started")
               if Tx>=TsAx && Tx<TsA[iTsA+1]
                   numofn[index]=(numofnA[iTsA+1]-numofnA[iTsA])/(TsA[iTsA+1]-TsAx)*(Tx-TsAx)+numofnA[iTsA]                   
                   #println("TsAx<Tx<TsAx+1 and numofn",numofn[index])
                   index=index+1
                   break
                end
                #println("iTsA between ended")
            end              
            #println("++++++++++++++++++++")
        end
        sigmae=0.0
        sigmah=0.0
        seebecke=0.0
        seebeckh=0.0
        kee=0.0
        keh=0.0
        #println("#########################")
        #println("T=$Tx")
        nx=numofn[previous]
        #println("nx ",nx," numofn[previous] ",numofn[previous] )
        bndstTx.var[1]=Tx    
        tauTOTTx.variables[2]=Tx       
        #println("n=$nx")
        xsx=xs[1]
        tauTOTTx.variables[7]=nx*1e6
        tauTOTTx.variables[9]=numofnn[1]*1e6           
        bndstTx.var[2]=xsx
        #println("x=$xsx")
        tauTOTTx.variables[4]=xsx
        types.updatebnstTx(bndst)
        Efs[i]=Fermilevel_n(nx,bndst,Tx,xmax)                
        for (l,band) in enumerate(bndst.bands)                    
            if l==1                    
                if band.effMass>0
                    #println("$Efs[i,j,k] ,$Tx")
                    presigmae=sigma(tauTOTTx,band,Efs[i],Tx)                            
                    preseebecke=seebeck_Nominator(tauTOTTx,band,Efs[i],Tx) 
                    prekee=keint(tauTOTTx,band,Efs[i],Tx) 
                    sigmae=isequal(presigmae,NaN) ? 0.0 : presigmae
                    seebecke=isequal(preseebecke ,NaN) ? 0.0 : preseebecke 
                    kee=isequal(prekee ,NaN) ? 0.0 : prekee 
                else
                    presigmah=sigma(tauTOTTx,band,Efs[i],Tx) 
                    preseebeckh=seebeck_Nominator(tauTOTTx,band,Efs[i],Tx)
                    prekeh=keint(tauTOTTx,band,Efs[i],Tx)
                    sigmah=isequal(presigmah,NaN) ? 0.0 : presigmah
                    seebeckh=isequal(preseebeckh ,NaN) ? 0.0 : preseebeckh 
                    keh=isequal(prekeh ,NaN) ? 0.0 : prekeh
                end                        
            else
                if band.effMass>0
                    temp0=sigma(tauTOTTx,band,Efs[i],Tx)     
                    #println("temp0=$temp0")                        
                    presigmae= bandpre==band ? presigmae : temp0#sigma(tauTOTTx,band,Efs[i,j,k],Tx)
                    #println("presigmae=$presigmae")
                    preseebecke= bandpre==band ? preseebecke : seebeck_Nominator(tauTOTTx,band,Efs[i],Tx)
                    prekee= bandpre==band ? prekee : keint(tauTOTTx,band,Efs[i],Tx)
                    sigmae+=isequal(presigmae,NaN) ? 0.0 : presigmae
                    seebecke+=isequal(preseebecke ,NaN) ? 0.0 : preseebecke 
                    kee+=isequal(prekee ,NaN) ? 0.0 : prekee 
                    #println("sigmaein=$sigmae")
                else
                    presigmah= bandpre==band ? presigmah : sigma(tauTOTTx,band,Efs[i],Tx)
                    preseebeckh= bandpre==band ? preseebeckh : seebeck_Nominator(tauTOTTx,band,Efs[i],Tx) 
                    prekeh= bandpre==band ? prekeh : keint(tauTOTTx,band,Efs[i],Tx)
                    sigmah+=isequal(presigmah,NaN) ? 0.0 : presigmah
                    seebeckh+=isequal(preseebeckh ,NaN) ? 0.0 : preseebeckh 
                    keh+=isequal(prekeh ,NaN) ? 0.0 : prekeh
                end
                    bandpre=band
            end
                #println("l = $l , sigmae= $sigmae ")
                #println("e l= ",l,"  ",seebecke)
                #println("h l= ",l,"  ",preseebeckh/presigmah)
        end       
        sigmaes[i]=sigmae
        sigmahs[i]=sigmah
        sigmas[i]=sigmae+sigmah
        seebecks[i]=(seebecke+seebeckh)/sigmas[i]            
        seebeckes[i]=(seebecke)/sigmaes[i]            
        seebeckhs[i]=(seebeckh)/sigmahs[i]
        kbis[i]=(sigmae*sigmah)/(sigmae+sigmah)*((seebecke/sigmae-seebeckh/sigmah)^2)*Tx                
        kees[i]=kee/q/q/Tx-((seebecke/sigmae)^2)*sigmae*Tx
        kehs[i]=keh/q/q/Tx-((seebeckh/sigmah)^2)*sigmah*Tx
        kes[i]=kees[i]+kehs[i]
        sigmas[i]=isequal(sigmas[i],NaN) ? 0.0 : sigmas[i]
        seebecks[i]=isequal(seebecks[i],NaN) ? 0.0 : seebecks[i]
        kbis[i]=isequal(kbis[i],NaN) ? 0.0 : kbis[i]
        kees[i]=isequal(kees[i],NaN) ? 0.0 : kees[i]
        kehs[i]=isequal(kehs[i],NaN) ? 0.0 : kehs[i]
        kes[i]=isequal(kes[i],NaN) ? 0.0 : kes[i]            
    end   
    return (numofn,Efs,sigmas,seebecks,kes,kees,kehs,kbis,sigmaes,sigmahs,seebeckes,seebeckhs)
    #return numofn
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
                    #println("sigma ",sigmae)
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
    defnumofnn=Array{Float64}(length(numofn))
    for (inn,) in enumerate(defnumofnn)
       inn=1e17 
    end
   return  electronicPropsiso(bndst,xs,Ts,numofn,defnumofnn,tauTOTTx,xmax)
end
function electronicPropsiso(bndst,xs,Ts,numofn,numofnn,tauTOTTx,xmax)
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
    sigmaes=Array{Float64}(length(Ts),length(numofn),length(xs))
    sigmahs=Array{Float64}(length(Ts),length(numofn),length(xs))
    seebeckes=Array{Float64}(length(Ts),length(numofn),length(xs))
    seebeckhs=Array{Float64}(length(Ts),length(numofn),length(xs))            
    kees=Array{Float64}(length(Ts),length(numofn),length(xs))
    kehs=Array{Float64}(length(Ts),length(numofn),length(xs))
    kes=Array{Float64}(length(Ts),length(numofn),length(xs))
    kbis=Array{Float64}(length(Ts),length(numofn),length(xs))
    Efs=Array{Float64}(length(Ts),length(numofn),length(xs))     
    for (i,Tx) in enumerate(Ts)
        #println("T=$Tx")
        bndstTx.var[1]=Tx    
        tauTOTTx.variables[2]=Tx 
        for (j,nx) in enumerate(numofn)
            #println("n=$nx")
            k=j
            xsx=xs[j]
            tauTOTTx.variables[7]=nx*1e6
            tauTOTTx.variables[9]=numofnn[j]*1e6
            #for (k,xsx) in enumerate(xs)                
                sigmae=0.0
                sigmah=0.0
                seebecke=0.0
                seebeckh=0.0
                kee=0.0
                keh=0.0
                bndstTx.var[2]=xsx
                #println("x=$xsx")
                tauTOTTx.variables[4]=xsx
                types.updatebnstTx(bndst)
                Efs[i,j,k]=Fermilevel_n(nx,bndst,Tx,xmax)                
                for (l,band) in enumerate(bndst.bands)                    
                    if l==1                    
                        if band.effMass>0
                   #         println("$Efs[i,j,k] ,$Tx")
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
                            temp0=sigma(tauTOTTx,band,Efs[i,j,k],Tx)     
                            #println("temp0=$temp0")                        
                            presigmae= bandpre==band ? presigmae : temp0#sigma(tauTOTTx,band,Efs[i,j,k],Tx)
                            #println("presigmae=$presigmae")
                            preseebecke= bandpre==band ? preseebecke : seebeck_Nominator(tauTOTTx,band,Efs[i,j,k],Tx)
                            prekee= bandpre==band ? prekee : keint(tauTOTTx,band,Efs[i,j,k],Tx)
                            sigmae+=isequal(presigmae,NaN) ? 0.0 : presigmae
                            seebecke+=isequal(preseebecke ,NaN) ? 0.0 : preseebecke 
                            kee+=isequal(prekee ,NaN) ? 0.0 : prekee 
                            #println("sigmaein=$sigmae")
                        else
                            presigmah= bandpre==band ? presigmah : sigma(tauTOTTx,band,Efs[i,j,k],Tx)
                            preseebeckh= bandpre==band ? preseebeckh : seebeck_Nominator(tauTOTTx,band,Efs[i,j,k],Tx) 
                            prekeh= bandpre==band ? prekeh : keint(tauTOTTx,band,Efs[i,j,k],Tx)
                            sigmah+=isequal(presigmah,NaN) ? 0.0 : presigmah
                            seebeckh+=isequal(preseebeckh ,NaN) ? 0.0 : preseebeckh 
                            keh+=isequal(prekeh ,NaN) ? 0.0 : prekeh
                        end
                        bandpre=band
                    end
                #println("l = $l , sigmae= $sigmae ")
               # println("e l= ",l,"  ",seebecke)
                #println("h l= ",l,"  ",preseebeckh/presigmah)
                end
                sigmaes[i,j,k]=sigmae
                sigmahs[i,j,k]=sigmah
                sigmas[i,j,k]=sigmae+sigmah
                seebecks[i,j,k]=(seebecke+seebeckh)/sigmas[i,j,k]            
                seebeckes[i,j,k]=(seebecke)/sigmaes[i,j,k]            
                seebeckhs[i,j,k]=(seebeckh)/sigmahs[i,j,k]
                kbis[i,j,k]=(sigmae*sigmah)/(sigmae+sigmah)*((seebecke/sigmae-seebeckh/sigmah)^2)*Tx                
                kees[i,j,k]=kee/q/q/Tx-((seebecke/sigmae)^2)*sigmae*Tx
                kehs[i,j,k]=keh/q/q/Tx-((seebeckh/sigmah)^2)*sigmah*Tx
                kes[i,j,k]=kees[i,j,k]+kehs[i,j,k]
                sigmas[i,j,k]=isequal(sigmas[i,j,k],NaN) ? 0.0 : sigmas[i,j,k]
                seebecks[i,j,k]=isequal(seebecks[i,j,k],NaN) ? 0.0 : seebecks[i,j,k]
                kbis[i,j,k]=isequal(kbis[i,j,k],NaN) ? 0.0 : kbis[i,j,k]
                kees[i,j,k]=isequal(kees[i,j,k],NaN) ? 0.0 : kees[i,j,k]
                kehs[i,j,k]=isequal(kehs[i,j,k],NaN) ? 0.0 : kehs[i,j,k]
                kes[i,j,k]=isequal(kes[i,j,k],NaN) ? 0.0 : kes[i,j,k]
            #end
        end
   end
    return (Efs,sigmas,seebecks,kes,kees,kehs,kbis,sigmaes,sigmahs,seebeckes,seebeckhs)
end

function thermalPropsiso(Ts,tauPHL::Array{tau_phonon_B,1},tauPHTx::Array{tau_phonon_B,1},tauPHTy::Array{tau_phonon_B,1},vsos::Array{Float64,1})
    
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
function differentialcond(band,xs,Ts,numofn,numofnn,tauTOTTx,xmax,Efs)
    bandpre=band
    Ex=collect(0.0:0.1:20.0)
    sigmaDs=Array{Float64}(length(Ts),length(numofn),length(xs),length(Ex))    
    #Efs=Array{Float64}(length(Ts),length(numofn),length(xs))     
    for (i,Tx) in enumerate(Ts)
        #println(rintln("T=$Tx")
        band.var[1]=Tx    
        tauTOTTx.variables[2]=Tx 
        for (j,nx) in enumerate(numofn)
            #println("n=$nx")
            k=j
            xsx=xs[j]
            tauTOTTx.variables[7]=nx*1e6
            tauTOTTx.variables[9]=numofnn[j]*1e6
            #for (k,xsx) in enumerate(xs)                
            sigmae=0.0                
            band.var[2]=xsx
            #println("x=$xsx")
            tauTOTTx.variables[4]=xsx
            #types.updatebnstTx(bndst)
            types.bandTxupdate(band)   
            for (Exi,Exx) in  enumerate(Ex)
                #Efs[i,j,k]=Fermilevel_n(nx,bndst,Tx,xmax)
                sigmaDs=sigmaD(tauTOTTx,band,Exx,Efs[i,j,k],Tx)   
                #sigmae=isequal(presigmae,NaN) ? 0.0 : presigmae              
            end                
        end
   end
    return (Ex,sigmaDs)
end
function kltSA(gammaSA,GM,Tts,MSiGecgs,thetaD,omegaD,beta,delta,Eep,mds,Efsr)
    klt_SA_A=Array{Float64,}(length(Ts))
    for Tti in 1:length(Tts)
        klt_SA_A[Tti]=klSA(gammaSA,GM,Tts[Tti],MSiGecgs,thetaD,omegaD,beta,delta,Eep,mds,Efsr[Tti,1,1]/kBe/Tts[Tti])   
    end
    return klt_SA_A
end


