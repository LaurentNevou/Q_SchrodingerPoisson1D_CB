%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Computes ISB dipoles %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Take care! Some people use meff inside the oscillator strenght f
% Actually, meff has sens in an infinite QW because there is a single mass value
% but NOT in multi-QW structure with various materials
% https://www.nextnano.com/nextnano3/tutorial/1Dtutorial_IntrabandTransitions.htm
% https://www.nextnano.com/nextnano3/tutorial/1Dtutorial_InGaAs_MQWs.htm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z_dipole_c = zeros(length(Ec),length(Ec));
f_dipole_c = zeros(length(Ec),length(Ec));

for i=1:length(Ec)
  for j=1:length(Ec)
    if j>i
      z_dipole_c(i,j) = abs(  trapz( z , psic(:,i).*z'.*psic(:,j) )  );
      f_dipole_c(i,j) = 2*m0/hbar^2 * ( Ec(j)-Ec(i) )* e * z_dipole_c(i,j)^2 ;
    end
  end
end

z_dipole_c = z_dipole_c + z_dipole_c.';
f_dipole_c = f_dipole_c + f_dipole_c.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%