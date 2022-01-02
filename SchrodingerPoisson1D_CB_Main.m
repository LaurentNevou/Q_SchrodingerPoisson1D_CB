%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% last update 25December2021, lne %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This program solves the Schrodinger-Poisson equation in the conduction band 
% for any heterostructures.
% The program solves the Schrodinger equation with m(E,z) using the Scaning/Shooting
% method. The non-parabolicity is implemented via the alpha parameter, alpha=1/Egap;
% meff(E)=meff(0)*(1+alpha*E); 
% It can be easily modified by tuning the alpha in each layer => line 109
% A strain model is included. It basically shifts the conduction band edge
% The strain is mainly interesting for InGaAs/GaAs heterostructures
% The non-parabolicity is also included into the density of states for the Poisson solver.
% It follows the book of Paul Harrison and actually makes the 2d density of states not constant
% -> Additionnal material can be added in the "materialDB_ZB.csv" file
% -> II-VI and cubic nitride material parameters are available but should
% be grabt in the "Library.m" file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If the code doesn t converge:
% -> decrease the doping
% -> increase the resolution dz
% -> increase the temperature (T=0K is very bad while T=10K is already much better)
% -> increase the amount of loops, Nloops
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h    = 6.62606896E-34;              %% Planck constant [J.s]
hbar = h/(2*pi);
e    = 1.602176487E-19;             %% electron charge [C]
m0   = 9.10938188E-31;              %% electron mass [kg]
Epsi0= 8.854187817620E-12;          %% Vaccum dielectric constant [F/m]
kB   = 1.3806488E-23;               %% Boltzmann's constant [J/K]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Nloops = 5;                   % number of loops
StrainModel = 0;              % Activate Strain model
n      = 4;                   % number of solution asked per model
ScF    = 0.1;                 % scaling factor to plot the wave function [Without Dimension]
dz     = 1e-10;               % resolution of the grid [m]
T      = 300;                 % Temperature [Kelvin]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DisplayResults   = 1;         % Switch to print or not the ISB dipoles on the shell

plot_density     = 1;         % Activate the plot 0 or 1
plot_convergence = 0;         % Activate the plot 0 or 1
plot_field       = 0;         % Activate the plot 0 or 1
plot_Vbending    = 0;         % Activate the plot 0 or 1
plot_mass        = 0;         % Activate the plot 0 or 1
plot_ro          = 0;         % Activate the plot 0 or 1
plot_Epsi        = 0;         % Activate the plot 0 or 1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Library;                  % load material parameter DB from "materialDB_ZB.csv"
ExtractParameters;        % extract parameter from the Library
TernaryAlloy;             % compute the ternary alloy    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% import the layer structure file %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the material used from the "library"
% second column is the length of the layer in nm
% third column is the n doping volumique of that layer in 1e18 cm-3

% You have to put a resonable amount of doping! Otherwise, it will diverge 
% and you will have to damp it more by increasing the number of loops !

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

input_file;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% NOTHING TO CHANGE ANYMORE !!! %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Grabbing the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

zt    = M(:,end-1)*1e-9;      % conversion of the length from nm to meter
Dopt  = M(:,end)*1e18*1e6;    % n doping conversion from cm-3 to m-3

Egt   = M(:,idx_Eg6c) - (M(:,idx_alphaG)*T^2) ./ (T+M(:,idx_betaG));   %Eg = Eg0 - (a*T.^2)./(T + b);
CBOt  = Egt+M(:,idx_VBO);     % CBO form band gap difference and temperature
%Dsot = M(:,idx_Dso);        % Spin-Orbit shift band parameter
%EPt_K= M(:,idx_EP_K);       % EP Kane
Epsit = M(:,idx_Epsi);        % Epsilon; dielectric constant

mefft = M(:,idx_me);          % electron effective mass at band edge
alphat= 1./Egt;               % non-parabolicity parameter => can be adjusted!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Strain Model %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

at  = M(:,idx_a);           % lattice parameter
act = M(:,idx_ac);          % Conduction band strain offset parameter
avt = M(:,idx_av);          % Valence band strain offset parameter
bvt = M(:,idx_bv);          % Valence band strain offset parameter
c11t = M(:,idx_c11);        % strain parameter
c12t = M(:,idx_c12);        % strain parameter

a0   = substrate(idx_a);

if StrainModel == 1
  exxt =  (a0-at)/a0; % eyyt =  exxt;
  ezzt = -2*c12t./c11t.*exxt;
else
  exxt =  (a0-at)/a0 * 0; % eyyt =  exxt;
  ezzt = -2*c12t./c11t.*exxt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Discretisation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here, I descretize the grid z, the potential V0 and that values that are needed

z=0; V0=CBOt(1); Eg=Egt(1); 
meff=mefft(1); alpha=alphat(1);
Dop=Dopt(1); Epsi=Epsit(1);
ac=act(1); av=avt(1); bv=bvt(1); exx=exxt(1); ezz=ezzt(1);

for i=1:length(zt)
    t    = zt(i);
    zv   = (z(end)+dz): dz : (z(end)+dz)+t;
    z    = [z zv];
    V0   = [ V0     ones(size(zv)) * CBOt(i)  ];
    Eg   = [ Eg     ones(size(zv)) * Egt(i)   ];
    Dop  = [ Dop    ones(size(zv)) * Dopt(i)  ];
    Epsi = [ Epsi   ones(size(zv)) * Epsit(i) ];
    meff = [ meff   ones(size(zv)) * mefft(i) ];
    alpha= [ alpha  ones(size(zv)) * alphat(i)];
    ac   = [ ac     ones(size(zv)) * act(i)   ];
    av   = [ av     ones(size(zv)) * avt(i)   ];
    bv   = [ bv     ones(size(zv)) * bvt(i)   ];
    exx  = [ exx    ones(size(zv)) * exxt(i)  ];
    ezz  = [ ezz    ones(size(zv)) * ezzt(i)  ];
end

V0=V0-min(V0);             % Shift the band in order to get the bottom of the well at zero
V0=(F0*z)+V0;              % adding the electric field to the potential

Ltot=z(end)-z(1);

eyy = exx;
DCBO   = -abs(ac).*(exx+eyy+ezz) ;                      % shift of the CB due to strain

Ntott=Dopt.*zt;
Ntot=sum(Ntott);   % total number of charges

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dE1=1e-3;
dE2=1e-2;
    
if StrainModel == 0
    E1 = min(V0):dE1:min(V0)+0.2 ;
    E2 = E1(end):dE2:max(V0);
    En=sort([E1 E2]);
else
    E1 = min(V0+DCBO):dE1:min(V0+DCBO)+0.2 ;
    E2 = E1(end):dE2:max(V0+DCBO);
    En=sort([E1 E2]);
end

EEn     = repmat(En',   [1  length(z)]);
V0mat   = repmat(V0 ,   [length(En) 1]);
Egmat   = repmat(Eg ,   [length(En) 1]);
meffmat = repmat(meff , [length(En) 1]);
alphamat= repmat(alpha, [length(En) 1]);
melin   = meffmat.*(1+alphamat.*(EEn-V0mat)); % => linearized mass

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Starting of the Poisson s loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Vs=zeros(size(z)); Vsold=Vs;
ntot=0; nloop=1;
ErrVec=1; sumVtotVec=1;

while nloop<Nloops
    
    nloop
    x = 0.5;
    Vbending=Vs*x + Vsold*(1-x);
    Vtot=V0+Vbending;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% schrodinger solver %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    dE=0.002; precision=1e-7;
    Vtotmat = repmat(Vtot, [length(En) 1]);
    melin = meffmat.*(1+alphamat.*(EEn-Vtotmat)); % => linearized mass
    [Ec,psic] = Schrod_Nbands_shoot_f(z,Vtot+DCBO,melin,n,En,dE,precision);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Here, I re-define the energy grid in order optimize the meshing
    dE1=1e-4;
    dE2=1e-2;
    if StrainModel == 0
        E1 = Ec(1):dE1:Ec(1)+0.1 ;
        E2 = E1(end):dE2:max(Vtot);
        En=sort([E1 E2]);
    else
        E1 = Ec(1):dE1:Ec(1)+0.1 ;
        E2 = E1(end):dE2:max(Vtot+DCBO);
        En=sort([E1 E2]);
    end

    EEn     = repmat(En' , [1 length(z) ]);
    V0mat   = repmat(V0  , [length(En) 1]);
    Vtotmat = repmat(Vtot, [length(En) 1]);
    Egmat   = repmat(Eg  , [length(En) 1]);
    meffmat = repmat(meff ,[length(En) 1]);
    alphamat= repmat(alpha, [length(En) 1]);
    melin   = meffmat.*(1+alphamat.*(EEn-Vtotmat)); % => linearized mass
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Paul Harrisson
    % Quantum Wells, Wires and Dots.
    % 4th edition (2016),
    % chap 2 : "Solutions to Schrodinger's equation"
    % 2.42: "Two-dimensional systems" page 31
    
    meffro = meffmat.*(1+2*alphamat.*EEn); % => linearized mass
    ro=[];
    for i=1:length(Ec)
        ro( En>Ec(i),:,i) = e*meffro(En>Ec(i),:) * m0/(pi*hbar^2);
        ro( En<Ec(i),:,i) = 0;
    end
    
    % here is what make the code very slow
    if Ntot==0
        idx=find(min(CBOt)==CBOt);
        Ef=-Egt(idx(1))/2;
        NN=0*Ec';
        roEf=ro*0;
    else
        [Ef,NN,roEf]=find_Ef_f(z,Ec,psic,En,ro,Ntot,T);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ntot2 = repmat(NN,[length(z) 1]).*abs(psic).^2 ;
    
    ntot = sum(ntot2,2)' - Dop;  % remove the charge positives (ions)
        
    %%%%%%%%%%%%%%%%%%%%%%%%%% Electrical Field %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    F = e*cumtrapz(z,ntot)./(Epsi0*Epsi);
    MF  = trapz(z,F)/(z(end)-z(1));  % MF=mean(F) on a nonlinear grid z
    F = F-MF;
    %%%%%%%%%%%%%%%%%%%%%%%%%%% New Potential %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Vsold = Vs;
    Vs    = -cumsum(F)*dz;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%% Convergence %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    Err = abs(  1 - sumVtotVec(end)/sum(Vs)  );
    sumVtotVec(nloop) = sum(Vs);
    ErrVec = [ErrVec Err];

    nloop=nloop+1;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Display Results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

scale_PSI;
computesISBdipoles;

if DisplayResults == 1
  PrintResults;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure('position',[-3500 100 1000 700],'color','w');
figure('position',[10 100 1000 700],'color','w');
subplot(1,1,1,'fontsize',15)
hold on; grid on; box on;
col=colormap(jet);

if plot_density==1
    grid off
    pcolor(z*1e9,En,sum(ROEf,3)*1e-6 )
    set(gca,'color',col(1,:))
    shading flat
    hcb=colorbar;
    caxis([0 max(max(sum(ROEf,3)*1e-6))+1])
    title(hcb,'\fontsize{8}cm-3')
    
    if StrainModel == 0
        plot(z*1e9,V0,  'w--','linewidth',1)
        plot(z*1e9,Vtot,'w-' ,'linewidth',1)
    else
        plot(z*1e9,V0+DCBO,  'k--','linewidth',1)
        plot(z*1e9,Vtot+DCBO,'k-' ,'linewidth',1)
    end
elseif plot_density==0
    if StrainModel == 0
        plot(z*1e9,V0,  'b--','linewidth',1)
        plot(z*1e9,Vtot,'b-' ,'linewidth',1)
    else
        plot(z*1e9,V0+DCBO,  'k--','linewidth',1)
        plot(z*1e9,Vtot+DCBO,'k-' ,'linewidth',1)
    end
end

for i=1:length(Ec)
    plot(z*1e9,PSIc(:,i),'color','r','linewidth',1)
end

plot([z(1) z(end)]*1e9,[1 1]*Ef,'g','linewidth',1)
text(z(end)*1e9*0.95,Ef+0.01,'\color{green}Fermi')

xlabel('z (nm)')
ylabel('Energy (eV)')
xlim([0 z(end)*1e9])
if StrainModel == 1
  title(strcat('T=',num2str(T),'K ; dz=',num2str(dz*1e9),'nm; Ntot=',num2str(Ntot*1e-4,'%.1e'),'cm-2; with STRAIN'))
else
  title(strcat('T=',num2str(T),'K ; dz=',num2str(dz*1e9),'nm; Ntot=',num2str(Ntot*1e-4,'%.1e'),'cm-2; without STRAIN'))
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_convergence==1
    figure('color','w')
    semilogy(1:nloop,ErrVec,'bo-')
    hold on; grid on; box on;
    xlabel('Cycles')
    ylabel('Convergence (norm. units)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_field==1
    
    figure('color','w')
    hold on; grid on; box on;
    [AX,H1,H2]=plotyy(z*1e9,F*1e-2*1e-3,z*1e9,Dop*1e-18*1e-6);
        
    set(H1,'color','r')
    set(H2,'color','b')
    
    xlabel('z (nm)')
    ylabel(AX(1),'E- field (kV/cm)','color','red')
    ylabel(AX(2),'Doping (1e18 cm-3)','color','blue')
    
    set(AX(1),'ycolor','red')
    set(AX(2),'ycolor','blue')
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_Vbending==1

    figure('color','w')
    hold on; grid on; box on;

    [AX,H1,H2]=plotyy(z*1e9,Vbending,z*1e9,ntot*1e-18*1e-6);

    set(H1,'color','r')
    set(H2,'color','b')
    
    xlabel('z (nm)')
    ylabel(AX(1),'Vbending (eV)','color','red')
    ylabel(AX(2),'ntot (1e18 cm-3)','color','blue')
    
    set(AX(1),'ycolor','red')
    set(AX(2),'ycolor','blue')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_mass==1
    figure('color','w')
    
    subplot(1,1,1)
    hold on; grid on; box on;
    pcolor(z*1e9,En,melin)
    
    colormap(jet)
    hcb=colorbar;
    title(hcb,'\fontsize{8}meff')
    shading flat
    xlabel('z (nm)')
    ylabel('Energy (eV)')
    title('me-lin')
    xlim([0 z(end)*1e9])
    
%     figure
%     hold on;grid on;
%     idx=250;
%     plot(En,melin(:,idx),'r')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_ro==1
    figure('color','w')
    hold on; grid on; box on;
    surf(z*1e9,En,sum(ro,3))
    view(45,30)
    colormap(jet)
    hcb=colorbar;
    title(hcb,'\fontsize{8}ro')
    shading flat
    xlabel('z (nm)')
    ylabel('Energy (eV)')
    xlim([0 z(end)*1e9])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if plot_Epsi==1
    figure('color','w')
    hold on; grid on; box on;
    plot(z*1e9,Epsi,'.-')
    xlabel('z (nm)')
    ylabel('Epsilon')
    xlim([0 z(end)*1e9])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
