function [states]=generateBands(lattdepth,waist,NSTM,NPTS,numPeriods)
% generateBands is used to calculate full wavefunctions in a lattice. 
    % The potential is given by the lattice depth (in recoils) with a
    % gaussian envelope with some waist (in microns), solved over
    % numPeriods periods of a 1064 lattice. This function will solve for
    % the first NSTM states, and discretizes the mesh (potential and
    % wavefunctions) into NPTS points.
    loadconstants;

    WFsum1D_1064s=[];
    gaussianwaist=waist*1e-6; %in microns
    depth1064=lattdepth*Er;

%     NSTM=200; % number of eigenvalues to calculate
%     NPTS=10000; % number of points for wavefunction discretization

    xmaxdiagonalize=numPeriods*1.064E-6;
    x=0:(xmaxdiagonalize/NPTS):xmaxdiagonalize; % in meters
    wavelength1064=1064E-9;

    twomoverhbarsquared=2*m_Sr/(hbar^2);

    latt1064=depth1064.*sin(2*pi*x./wavelength1064).^2;
    
    gauss1064=-exp(-2*(x-mean(x)).^2/(gaussianwaist^2));
    V_1064=gauss1064.*latt1064;
    V_1064=V_1064-min(V_1064);
%     figure(1);
%     plot(x, gauss1064); hold on;
%     plot(x, V_1064);

    %figure(88);clf;subplot(211);
    %plot(x,latt1550,'r','LineWidth',2);hold on;
    %set(gcf,'Color','white');
    %plot(x,latt775,'b','LineWidth',2);
    %plot(x,latt775up,'r:','LineWidth',1);
    %plot(x,latt775dn,'b:','LineWidth',1);

    %subplot(212);
    %plot(x,latt1550+latt775,'m:','LineWidth',3); hold on;
    %plot(x,latt1550+latt775up,'c--','LineWidth',3);
    %plot(x,latt1550+latt775dn,'g--','LineWidth',3);

    % now find the eigenvalues

    %x=0:.1:10;
    %y=30*(latt1550+latt775up);
    potential=twomoverhbarsquared*V_1064;
    matty=[x;potential];
    dlmwrite('potdat3',matty);
    [ee, ev, V,x]=qm1d_fast_dmw('potdat3',NPTS,NSTM,xmaxdiagonalize);
    groundphi=-ev(:,1)';

    figure(2);clf;
    % get ee in recoil units
    ee=ee.*hbar^2./(2*m_Sr);
    ee=ee./Er;
    % average spacing of energy levels (for adjusting scale of ev)
    de0=(ee(NSTM,NSTM)-ee(1,1))/(NSTM-1);
    % plot potential
    Vplot=V./twomoverhbarsquared; Vplot=Vplot./Er;
    plot(x,Vplot,'k'); hold on;
    % plot eigenvectors
    de=0.35*sqrt(NPTS)*de0;
    for n=1:NSTM
        plot(x,ee(n,n)+de*(-1)*ev(:,n),'b'); hold on;
        states(n,:) = de*(-1)*ev(:,n); %+ee(n,n)
    end
    xlim([0 xmaxdiagonalize]); ylim([min(Vplot) max(Vplot)]); set(gcf,'Color','white');%ee(NSTM,NSTM)+de0

    % % now calculate U
    % dx=xmaxdiagonalize/NPTS;
    % groundphinorm=groundphi/sqrt(dx);
    % WFsum=(sum(groundphinorm.^4)*dx)^3;
    % WFsum1D_1064=sum(groundphinorm.^4)*dx;
    % U=(4*pi*a11*hbar^2/m_Sr)*WFsum; % see PRL 81, 3108 for this expression for U

    % % now make the longer versions of X, V, and groundphi so we can do the
    % % nearest-neighbor integral
    % zeropad=zeros(1,length(x(x<(wavelength1064/2))));
    % WF1=[groundphi zeropad];
    % WF2=[zeropad groundphi];
    % xpad=zeropad; for boo=1:length(xpad) xpad(boo)=x(boo); end 
    % xlonger=[x (xpad+max(x))];
    % Vpad=zeropad; for boo=1:length(Vpad) Vpad(boo)=V(boo); end
    % Vlonger=[V Vpad]; % note that this assumes that V contains an exact number of half-wavelengths.

    % figure(2);clf;subplot(211);
    % plot(x,V,'k','LineWidth',3); hold on;
    % plot(x,3*max(V)*groundphi,'b','LineWidth',2); 
    % plot(x+wavelength1064/2,3*max(V)*groundphi,'g','LineWidth',2);
    % set(gcf,'Color','white');
    % 
    % J=getJ(xlonger,WF1,WF2,Vlonger);
    