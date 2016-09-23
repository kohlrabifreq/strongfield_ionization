function [J] = getJ(X,WF1,WF2,V)

% this function is supposed to calculate the tunneling matrix element given two
% adjacent wavefunctions and the potential.  it's modeled after qm1d_fast_dmw.m

% by the time this function takes them in, all these vectors should be the
% same length.

loadconstants;
NPTS=length(X);
L=max(X);

j=1:NPTS; % indexes for main diagonal
h=L/(NPTS-1); % space step

main_diag=2/h^2+V(j);
sub_diag=-1/h^2*ones(1,NPTS-1);
d=[0 -1 1];
B=[main_diag' [0; sub_diag'] [sub_diag'; 0]];
AA=spdiags(B,d,NPTS,NPTS); % create a sparse tridiagonal matrix representing the hamiltonian.

J=(hbar^2./(2*m_Rb)).*WF1*AA*WF2';

