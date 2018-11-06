function uexac=solexac_aPVI(t)    
    uexac(1,:)=  2*exp(-t)+sin(t);
    uexac(2,:) =2*exp(-t)+cos(t) ;
end