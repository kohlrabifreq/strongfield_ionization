% this code calculates the interaction energy in a simple 3D 1064nm lattice
% at different lattice depths.  the calculation is done using a finite difference method,
% which means matlab has to diagonalize a very large tridiagonal matrix.
% the code then plots the results of that calculation along with the
% results of a harmonic oscillator approximation and the numerical results
% from the Jaksch et al paper (PRL 81, 3108)

clear;

loadconstants;

Us=[];
Js=[];
WFsum1D_1064s=[];
lattdepth=1:20;

for baba=lattdepth
    depth1064=baba*Er;
    diagonalizelatticeH_1064;
    Us=[Us U]; % this is my calculation of U
    Js=[Js J]; % this is my calc of J
    WFsum1D_1064s=[WFsum1D_1064s WFsum1D_1064];
end
figure(3); clf; subplot(121);
set(gcf,'Color','white');
set(gca,'FontSize',14);
plot(lattdepth,Us./h,'b','LineWidth',2); hold on;
U=sqrt(8/pi)*k*a11*Er*(lattdepth.^(3/4)); % this is the SHO approximation for U
plot(lattdepth,U./h,'r--','LineWidth',2);
% following commented out by RS as they refer to non-existent file
% cd /Users/dweld/Documents/physics/code/kett'erle matlab code'/
% a=dlmread('jakschcalcsU.txt');
% lattdepthjaksch=a(:,1);
% Ujaksch=a(:,2)*Er*a11/532E-9/h; % these data were read in from the graph in the Jaksch et al PRL
% plot(lattdepthjaksch,Ujaksch,'k:','LineWidth',2);
xlabel('Lattice Depth (E_R)');
ylabel('Onsite Interaction Energy (Hz)');
legend('BEC4 Calculation','DDL Approximation','Location','NW'); %RS removed Jaksch label

subplot(122);
set(gca,'FontSize',14);
semilogy(lattdepth(lattdepth>4),Js(lattdepth>4)/h,'b','LineWidth',2); hold on;
JapproxDDL=(4/sqrt(pi))*(Er).*(lattdepth.^0.75).*exp(-2.*sqrt(lattdepth)); % this is the DDL approximation for J
plot(lattdepth,JapproxDDL/h,'r--','LineWidth',2);

% The following were commented out by RS as they refer to a non-existent
% file
% a=dlmread('jakschcalcsJ.txt');
% lattdepthjaksch2=a(:,1);
% Jjaksch=a(:,2)*Er/h/1000; % these data were read in from the graph in the Jaksch et al PRL.  The 1000 is because I read them in using the other y axis.
% plot(lattdepthjaksch2,Jjaksch,'k:','LineWidth',2);
% JapproxAHDL=(pi^2/4)*lattdepth.*Er.*exp(-(pi^2/4).*sqrt(lattdepth));
% plot(lattdepth,JapproxAHDL/h,'m.-','LineWidth',2);
%JapproxHM=3*sqrt(lattdepth.*Er.*Er).*exp(-(pi^2/4).*sqrt(lattdepth));
%plot(lattdepth,JapproxHM/h,'g','LineWidth',1);

xlabel('Lattice Depth (E_R)');
ylabel('Tunneling Energy (Hz)');
legend('BEC4 Calculation','DDL Approximation','AHDL Approximation','Location','NE'); %RS removed Jaksch label

save('UandJdata_1064','Js','WFsum1D_1064s','lattdepth');