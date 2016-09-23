% This routine does the following:
% 1. Reads the potential from input file.
% 2. Do spline interpolation for NPTS points.
% 3. Discretizes the schrodinger equation over the specified set of points:
%     v(x)=2m*V(x)/hbar^2; ee=2m*E/hbar^2.
%     [    d^2       2m        ]           2m*E
%     [ - -----  + ------ V(x) ] phi(x) = ------  phi(x)
%     [   d x^2    hbar^2      ]          hbar^2
%
% 4. Finds eigenvalues and eigenvectors.
% 5. Plots em.
%
% Usage:
% [ee,ev] = qm1d_fast(pot_filename,NPTS,NSTM,L)
% pot_file - string that specifies filename with the potential
% NPTS - number of points for discretization of schrodinger equation
% NSTM - number of eigen values and eigen vector to find
% L - the length of the interval, starting from zero.
%
% Example: qm1d_fast('pot.dat',10000,5,10);
%
% See also: http://iffwww.iff.kfa-juelich.de/~ekoch/DFT/qm1d.html
function [ee,ev,V,x] = qm1d_fast_dmw(pot_filename,NPTS,NSTM,L)
pot_dat=load(pot_filename);
j=1:NPTS; % indexes for main diagonal
h=L/(NPTS-1); % space step

x=j*h;
V=spline(pot_dat(1,:),pot_dat(2,:),x);

main_diag=2/h^2+V(j);
sub_diag=-1/h^2*ones(1,NPTS-1);

d=[0 -1 1];
B=[main_diag' [0; sub_diag'] [sub_diag'; 0]];
AA=spdiags(B,d,NPTS,NPTS); % create a sparse tridiagonal matrix
%AA
%[ee,ev] = trideigs(main_diag,sub_diag,'I',1,NSTM);
[ev,ee]=eigs(AA,NSTM,'sm');


