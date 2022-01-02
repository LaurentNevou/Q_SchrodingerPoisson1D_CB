function[Ef,NN,roEf]=find_Ef_f(z,Ec,psic,E,ro,Ntot,T)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

e  = 1.602176487E-19;           %% electron charge [C]
kB = 1.3806488E-23;             %% Boltzmann's constant [J/K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EE   = repmat(E'   ,[1 length(z)]);
   
%%%%%%%%%%%%%%%%%%% Computes the Fermi level at any T %%%%%%%%%%%%%%%%%%%%%%%%%%
if T==0
   T=1e-10; 
end
 
Ef=Ec(1);
Fermi= 1./(1+exp((EE-Ef)/(kB*T/e))) ;
FFermi=repmat(Fermi,[1 1 length(Ec)]);
roEf = ro.*FFermi;
NNz = squeeze(trapz(E,roEf,1));
NN = trapz(z,NNz.*abs(psic).^2 ,1)  ;
NtotX=sum(NN);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
% Now, it will try to get as close as posible to the real Ef with an
% error of 0.1% by dichotomy
ddE=0.01; % eV
Ef1=Ef;
Ef2=Ef1+ddE;

while  abs(NtotX - Ntot)/Ntot > 0.001  % find the Fermi level at any temperature
    
    if NtotX > Ntot
        Ef  = Ef - abs(Ef1-Ef2)/2 ;  
        Ef1 = Ef ;
    else
        Ef  = Ef + abs(Ef1-Ef2)/2 ;
        Ef2 = Ef ; 
    end
    
    Fermi  = 1./(1+exp((EE-Ef)/(kB*T/e))) ;  % Fermi Dirac distribution function
    FFermi = repmat(Fermi,[1 1 length(Ec)]);
    roEf   = ro.*FFermi;
    NNz    = squeeze(trapz(E,roEf,1));
    NN     = trapz(z,NNz.*abs(psic).^2 ,1) ;
    NtotX  = sum(NN);

end

end