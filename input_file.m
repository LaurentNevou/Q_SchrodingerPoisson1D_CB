%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Layers Structure %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first column is the material used from the "library"
% second column is the length of the layer in nm = 1e-9m
% third column is the doping of the layer in 1e18 cm-3 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

substrate=InP;
F0=0;%6e6;%0;               % Electric field [Volt/meter]

M=[
AlInAs      5   0
AlInAs      1   7
AlInAs      5   0
InGaAs     6   0
AlInAs      5   0
InGaAs     7   0
AlInAs      5   0
AlInAs      1   5
AlInAs      5   0
];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% substrate=GaAs;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
% F0=0;%6e6;%0;               % Electric field [Volt/meter]
% M=[
% GaAs      10  0.5
% InGaAs20  15  0
% GaAs      10  0.5
% ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% substrate=GaAs;      % Important for the Strain model (Si, GaAs, InP, InAs, GaSb)
% F0=0;%16e6;%0;               % Electric field [Volt/meter]
% M=[
% AlGaAs40      5  2
% GaAs         10  0
% AlGaAs40      5  2
% ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%