% Universal constants
hbar = 1.054E-34; % hbar in m^2 Kg / s
h = 2*pi*hbar;
c = 2.9979E8; % speed of light in m /s
mu_o= 4*pi*1E-7; % magnetic constant
a_o=0.53E-10; %bohr radius
kB=1.38E-23; % boltzmann constant
mu_B=9.27E-24; % bohr magneton

% % Sr constants -- new data from S. Kokkelmans revises these numbers!
% m_Rb = 88*1.67*10^-27; % atomic mass in Kg
% % a11=100.44*a_o; % scattering lengths (old values)
% % a12=98.09*a_o;
% % a22=95.47*a_o;
% a11=-1*a_o; % scattering lengths (new from servaas kokkelmans)
% %a12=98.98*a_o;
% %a22=98.98*a_o;

%Li constants
m_Sr=84*1.67*10^-27; %atomic mass in Kg
a11=-27*a_o;

% Experiment constants
lambda = 1064E-9; % lattice laser wavelength in m
k=2*pi/lambda;
Er = h^2/(2*m_Sr*lambda^2);
w_trap = 2*pi*155; % geometric mean of three trap frequencies in cyc/sec(needs to be updated!)
w_latt=74E-6; % geometric mean of lattice 1/e^2 radii (needs to be updated)