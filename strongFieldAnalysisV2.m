function [outputBand]=strongFieldAnalysisV2(lattdepth,freq)
    %strongFieldAnalysis(lattdepth, freq)
    % This function applies a pulse to atoms in a lattice
    %   -lattdepth is depth of lattice in recoils
    %   -freq is the frequency of the driving pulse in kHz
    %   -atoms are loaded into the ground state of the first band
    %   -the function shows a graph of the percentage of the wavefunction
    %   that remains in the first band
    
    %% Simulation Parameters and Things.
%     alphaval=str2double(alphaval);
%     Omega=str2double(Omega);
%     
%     a=8;
%     b=a*alphaval; 
%     lattdepth=25; %in recoils

    [hbar,h,c,mu_o,a_o,kB,mu_B,m_Sr,lambda,k,Er,w_trap]=loadconstantsfn;

    Omega=2*pi*freq*1e3;
    waist=50; %in microns
    gaussianwaist=waist*1e-6; %in m
    depth1064=lattdepth*Er;
    %w0 = sqrt((2*lattdepth*Er*k^2)/m_Sr);

    NSTM=20; % number of eigenvalues to calculate
    numPeriods=5;

    tstart=now;   
    %Position Mesh
    xmax=numPeriods*1.064E-6;
    NPTS=1000; % number of points in mesh
    dx = xmax/NPTS; % spacing between sites
    x=0:dx:xmax; % in meters
    wavelength1064=1064E-9;
    phi = 2*pi*x./wavelength1064;
    dphi = 2*pi*dx/wavelength1064;
    positionMesh=linspace(0,xmax,NPTS);%Make it with an extra
    %positionMesh=positionMesh(1:end-1);%Remove the boundary PBC
    x=x(1:end-1);

    %Time and Drive Mesh
    numDrive=100;
    maxTheta = pi;
    dTheta = maxTheta/numDrive;
    theta1 = 0:dTheta*2:maxTheta;
    theta2 = maxTheta-dTheta:-dTheta*2:0;
    thetaMesh = [theta1 theta2];
    totalTimesteps = length(thetaMesh); 
    
    dtau=dTheta./Omega; 
    dtau=dtau*hbar/(2*m_Sr);
    %Simulation Parameter
    kappa=2/(sqrt(lattdepth)*dx^2);
    simParam=kappa*dtau;
    %disp(['SimParam=' num2str(simParam)]);

    %% Create time evolution operators
    %Now make the single period time evolution operator.
    Tmatrix=gallery('tridiag',NPTS,-1,2,-1);%Make a sparse tridiagonal matrix    
 %  Tmatrix=Tmatrix+sparse([1 NPTS],[NPTS 1],[-1 -1],NPTS,NPTS);%add PBC
    %Tmatrix=kappa*Tmatrix;
    Tmatrix=Tmatrix./(dx^2);
    Imatrix=speye(NPTS);
    UDrivePeriod=Imatrix;
   
    %% Make Static Hamiltonian
    staticVmatrix=sparse(diag(staticpotential(depth1064,x,wavelength1064,gaussianwaist)));
    staticHmatrix=Tmatrix+staticVmatrix;
 
   %% Make UCrank matrix
    for ii=1:length(thetaMesh)        
        theta=thetaMesh(ii);
        Vmatrix=sparse(diag(relpotential(depth1064,x,wavelength1064,gaussianwaist,theta))); 
        Vmatrix=Vmatrix*2*m_Sr./(hbar^2);
        Hmatrix=Tmatrix+Vmatrix;
        UCrankNow=(Imatrix-1i*Hmatrix*(dtau/2))\(Imatrix+1i*Hmatrix*(dtau/2));
        UCrank(:,:,ii)=full(UCrankNow);
%         plot(x,relpotential(depth1064,x,wavelength1064,gaussianwaist,theta));
%         drawnow
    end
    
    %keyboard
    
    %% Initialize wavefuctions and time Evolve
%     init_WF_Up=makeGaussianComb(1);
%     init_WF_Down=makeGaussianComb(0);
    states=generateBands(lattdepth, waist, NSTM, NPTS, numPeriods);
    init_WF = states(1,:);
    
    norm = sum(init_WF.*conj(init_WF));
    init_WF = sqrt(1/norm)*init_WF.';
    rhoInit = init_WF.*conj(init_WF);
    states = states*sqrt(1/norm);
    figure(1);clf;
    plot(x,rhoInit,'LineWidth',3);hold on;
    
    %% Create Overlap Vectors
    rhos = [];
    for jj=1:NSTM
        rhos(jj,:) = states(jj,:).*conj(states(jj,:));
    end
    
    %% Perform End Analysis
    projections = [];
    outputBand=endAnalysis(init_WF);
    tend=now;
    
    disp('************');
    disp(['Elapsed Time : ' num2str(24*60*60*(tend-tstart))]);
    disp(outputBand);
    disp('************');   
    filename=['data_lattdepth' num2str(lattdepth) 'recoils_freq' num2str(freq) 'kHz.mat'];
    save(filename,'outputBand');
        
    %% Helper Functions.   

    function output=endAnalysis(psi)
        psi_i=psi;
        init_proj_data=calcOverlap2(psi_i);
        
       numCycles=1;
       Edata=[];%zeros(numCycles*numDrive,1);
       proj_data=[];%zeros(numCycles*numDrive,1);
       psi_data=[];%zeros(length(positionMesh),numCycles*numDrive);
       index=1;
       figure(4);clf;
       for ll=1:numCycles
          for jjj=1:totalTimesteps
              psi=UCrank(:,:,jjj)*psi;
              rho_now = psi.*conj(psi);
              plot(x,rho_now);
              drawnow
              %pause(.1);
              index=index+1;
          end
       end
       Edata=staticEnergy(psi);
       proj_data=calcOverlap2(psi);
       band_data=determineBandProb(psi);
       psi_data=psi;
       rhoFinal=psi.*conj(psi);
       output=struct;
       output.LatticeDepth=lattdepth;
       output.FreqkHz=freq;
       output.NumPerDriveCycle=numDrive;
       output.NumDriveCycles=numCycles;
       output.FinalWF = psi_data;
       output.StateProjections = proj_data;
       output.FirstBandPercent = band_data;
       output.FinalEnergy = Edata;
       figure(1);
       plot(x,rhoFinal,'r', 'LineWidth',1);
       
       %keyboard
       
%        output.Final4CycleEnergy=mean(Edata(end-numDrive*4:end));
%        output.Final4CycleDownStab=mean(omin_data(end-numDrive*4:end));
%        output.Final4CycleUpStab=mean(omax_data(end-numDrive*4:end));
%        output.Final4CycleStability=output.Final4CycleUpStab-output.Final4CycleUpStab;
%        output.WF_data=psi_data(:,1:50:end);       
    end

    function theans=calcOverlap(psi)
        rho=psi.*conj(psi);
        for kk=1:NSTM
            iterState = rhos(kk,:);
            iterState = iterState.';
            projections(kk) = sum(rho.*iterState);
            theans=projections;
        end
    end

    function theans=calcOverlap2(psi)
        theprojections=zeros(NSTM,1);
        for kk=1:NSTM
            thisWF=states(kk,:);
            thisnorm=sum(sum(thisWF.*conj(thisWF)));
            thisWF=thisWF'./sqrt(thisnorm);
            theprojections(kk) = sum(psi.*conj(thisWF));
        end
        theans=abs(theprojections);
    end

    function band1prob=determineBandProb(psi)
        projections=calcOverlap2(psi);
        if NSTM >10
            band1prob=sum(projections(1:10));
            band2prob=sum(projections(11:NSTM));
        else
            band1prob=sum(projections(1:NSTM));
        end
    end
        

    function E=staticEnergy(psi)
       E=abs(ctranspose(psi)*staticHmatrix*psi); 
    end
    
    function [val]=staticpotential(depth1064, x, wavelength1064, gaussianwaist)
       latt1064 = depth1064.*sin(2*pi*x/wavelength1064).^2;
       gauss1064=-exp(-2*(x-mean(x)).^2/(gaussianwaist^2));
       V1064=latt1064.*gauss1064;
       val=V1064-min(V1064);
    end

    function [val]=relpotential(depth1064, x, wavelength1064, gaussianwaist,phase)
       latt1064 = depth1064.*sin((2*pi*x/wavelength1064)+phase).^2;
       gauss1064=-exp(-2*(x-mean(x)).^2/(gaussianwaist^2));
       V1064=latt1064.*gauss1064;
       val=V1064-min(V1064);
    end

%     function [psi_i] = makeGaussian(center)
%     [psi_i]=arrayfun(@(phi) exp(-(sqrt(2*a)/8)*(phi-2*pi)^2),positionMesh);
%     normal=sqrt(sum(abs(psi_i).*abs(psi_i)));
%     
%     [psi_i]=arrayfun(@(phi) exp(-(sqrt(2*a)/8)*(phi-center)^2),positionMesh);
% 
%     psi_i=psi_i/normal;
%     psi_i=transpose(psi_i);
%     end
% 
%     function [psi]=makeGaussianComb(state)         
%        centers=-1:ceil(rangePhi/pi)+1;%Goes from -1 to N
%        centers=2*pi*centers;%Now gos from -2pi to N*2pi
%        centers=centers+pi*state;%if state=1 now goes from -pi to N*2pi-pi
%         %if state=1 inverted
%         %if state=0 nonverted
%        psi=zeros(length(positionMesh),1);
%        for index=1:length(centers);
%            center=centers(index);
%            thisState=makeGaussian(center);
%            psi=psi+thisState;
%        end
%        
%         normal=sqrt(sum(abs(psi).*abs(psi)));
%         psi=psi/normal;        
%         
%     end
    
end

