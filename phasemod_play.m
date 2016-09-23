    loadconstants;

    WFsum1D_1064s=[];
    lattdepth = 25; % in recoils
    waist = 50; %in microns
    gaussianwaist=waist*1e-6; %in meters
    depth1064=lattdepth*Er;
    numPeriods = 25;

    NSTM=200; % number of eigenvalues to calculate
    NPTS=10000; % number of points for wavefunction discretization

    xmaxdiagonalize=numPeriods*1.064E-6;
    x=0:(xmaxdiagonalize/NPTS):xmaxdiagonalize; % in meters
    wavelength1064=1064E-9;

    maxTheta = pi;
    numTheta = 100;
    dTheta = maxTheta/numTheta;
    theta1 = 0:dTheta:maxTheta;
    theta2 = maxTheta-dTheta:-dTheta:0;
    theta = [theta1 theta2];
    totalTimesteps = length(theta);
   
    %F(numTheta) = struct('cdata', [], 'colormap', []);
    
    for ii=1:totalTimesteps
        latt1064=depth1064.*sin((2*pi*x./wavelength1064)+theta(ii)).^2;
        gauss1064=-exp(-2*(x-mean(x)).^2/(gaussianwaist^2));
        V_1064=gauss1064.*latt1064;
        V_1064=V_1064-min(V_1064);
        plot(x, V_1064);
        drawnow
    end
    
        