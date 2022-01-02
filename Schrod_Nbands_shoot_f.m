function[E,psi]=Schrod_Nbands_shoot_f(z,V0,me,n,Evec,dE,precision)

method=2;         % method "2" is much more acurate than "1" but take also more time...
e=min(V0);
Emax=max(V0)+0.1;

C=0; N=0; E=[];
psi=z*0+1; psi_old=psi;


while e<Emax && length(E)<n
    
    C = C+1;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In this block, it is scanning to find the mass at the right Energy, m(E)

    idx=[0 0];
    epsi=1;
    while length(idx)>1
        idx=find( abs( (Evec-e )) < epsi ) ;
      if isempty(idx)
          idx=IDX(1);
        else
          IDX=idx;
          epsi=epsi/2;
        end
    end
    Mass=me(idx,:);
%    [ e Mass(round(end/2)) ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    psi_old = psi;
    psi=Schroed1D_Euler_Eval(z,V0,Mass,e,method);
    
    if (sign( psi(end) ) ~= sign( psi_old(end) ) ) && C>1 % here, I catch a quantum state because the last point of psi change sign
    
        N=N+1; de=dE;
        
        while abs(de)>precision
            if sign( psi(end) ) ~= sign( psi_old(end) )
                 de = -de/2;
            end
            e=e+de;
            psi_old=psi;
            psi=Schroed1D_Euler_Eval(z,V0,Mass,e,method);
        end

        E(N,:)=e; PSI(:,N)= psi;
        C=0;
    end	
    
    e=e+dE;
end

psi=PSI;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Normalization of the Wavefunction %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:N
    z=z';
    psi(:,i)=psi(:,i)/sqrt(trapz(z,abs(psi(:,i)).^2));  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[psi]=Schroed1D_Euler_Eval(z,V0,Mass,ee,method)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=6.62606896E-34;               %% Planck constant [J.s]
hbar=h/(2*pi);
e=1.602176487E-19;              %% electron charge [C]
m0=9.10938188E-31;              %% electron mass [kg]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dz=z(2)-z(1);
dz=[ dz diff(z)];
y = [0 1];

if method==1
    
    for i = 1:length(z)                                  %% number of z steps to take
        dy(2) = -2*e/(hbar^2) * ( ee  - V0(i) ) * y(1) ; %% Equation for dv/dz
        dy(1) =  Mass(i)*m0*y(2);                        %% Equation for dx/dz
        y = y + dz(i)*dy ;                               %% integrate both equations with Euler
        psi(i)= y(1);
    end
    
elseif method==2
    
    for i = 1:length(z)                                  %% number of z steps to take
    
        dy(2) = -2*e/(hbar^2) * ( ee  - V0(i) ) * y(1) ; %% Equation for dv/dz
        dy(1) =  Mass(i)*m0*y(2);                        %% Equation for dx/dz
        
        K = y + 0.5*dz(i)*dy;
        dK(2) = -2*e/(hbar^2) * ( ee  - V0(i) ) * K(1) ; %% Equation for dv/dz
        dK(1) =  Mass(i)*m0*K(2);                        %% Equation for dx/dz
        
        y =  y + dz(i)*dK;
        psi(i)= y(1);
    end
end

end