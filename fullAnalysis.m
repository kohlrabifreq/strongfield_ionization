
function fullAnalysis(alphaval,Omega)
    %singleStability(alphaval, Omega)
    %This function evaluates the stability of a point at a given alpha and
    %omega.
    %   alphaval: value of the alpha parameter
    %   Omega: value of the Omega parameter
    %   overlapMax: Projection of final state onto time average energy
    %   maximum
    %   overlapMin: Projection of final state onto time averged energy
    %   minimum
    
    %% Simulation Parameters and Things.
     alphaval=str2double(alphaval);
     Omega=str2double(Omega);
    
    a=8;
    b=a*alphaval; 
    
    tstart=now;   
    %Position Mesh
    numPhi=300;%Number of points in mesh
    rangePhi=6*pi;%Range
    positionMesh=linspace(0,rangePhi,numPhi+1);%Make it with an extra
    positionMesh=positionMesh(1:end-1);%Remove the boundary PBC
    dphi=rangePhi/numPhi;%Spacing between sites.

    %Time and Drive Mesh
    numDrive=500;
    thetaMesh=linspace(0,2*pi,numDrive+1);
    thetaMesh=thetaMesh(1:end-1);
    thetaMesh=thetaMesh+thetaMesh/numDrive;    
    UCrank=zeros(numPhi,numPhi,numDrive); 
    dtheta=2*pi./numDrive; 
    
    dtau=dtheta./Omega; 
    
    rangeTau=25*pi; 

    %Simumlation Parameter
    kappa=2/(sqrt(a)*dphi^2);
    simParam=kappa*dtau;
    %disp(['SimParam=' num2str(simParam)]);

    %% Create time evolution operators
    %Now make the single period time evolution operator.
    Tmatrix=gallery('tridiag',numPhi,-1,2,-1);%Make a sparse tridiagonal matrix    
    Tmatrix=Tmatrix+sparse([1 numPhi],[numPhi 1],[-1 -1],numPhi,numPhi);%add PBC
    Tmatrix=kappa*Tmatrix;
    Imatrix=speye(numPhi);  
    UDrivePeriod=Imatrix;
    
    %% Make Static Hamiltonian
    staticVmatrix=sparse(diag(arrayfun(@(phi) staticpotential(phi), positionMesh)));      
    staticHmatrix=Tmatrix+staticVmatrix;
    
   %% Make UCrank matrix
    for ii=1:length(thetaMesh)        
        theta=thetaMesh(ii);
        Vmatrix=sparse(diag(arrayfun(@(phi) relpotential(phi,theta), positionMesh)));      
        Hmatrix=Tmatrix+Vmatrix;
        UCrankNow=(Imatrix-1i*Hmatrix*(dtau/2))\(Imatrix+1i*Hmatrix*(dtau/2));
        UCrank(:,:,ii)=UCrankNow;
    end
    keyboard;
    %% Initialize wavefuctions and time Evolve
    init_WF_Up=makeGaussianComb(1);
    init_WF_Down=makeGaussianComb(0);

    
    %% Create Overlap Vectors
    combWF_Up=makeGaussianComb(1);
    combWF_Down=makeGaussianComb(0);    
    rho_min=abs(combWF_Down).*abs(combWF_Down);
    rho_max=abs(combWF_Up).*abs(combWF_Up);

    %% Perform End Analysis
    output_Up=endAnalysis(init_WF_Up);
    output_Down=endAnalysis(init_WF_Down);
    tend=now;
    
    disp('************');
    disp(['Elapsed Time : ' num2str(24*60*60*(tend-tstart))]);
    disp(output_Up);
    disp(output_Down);
    disp('************');   
    filename=['data_a' num2str(alphaval) '_b' num2str(Omega) '.mat'];
    save(filename,'output_Up','output_Down');
        
    %% Helper Functions.   

    function output=endAnalysis(psi)
       numCycles=300;
       Edata=zeros(numCycles*numDrive,1);
       omax_data=zeros(numCycles*numDrive,1);
       omin_data=zeros(numCycles*numDrive,1);
       psi_data=zeros(length(positionMesh),numCycles*numDrive);
       index=1;
       for ll=1:numCycles
          for jj=1:numDrive
              Edata(index)=staticEnergy(psi);
              omax_data(index)=maxOverlap(psi);
              omin_data(index)=minOverlap(psi);           
              psi_data(:,index)=psi;              
              psi=UCrank(:,:,jj)*psi;
              index=index+1;
          end
       end       
       output=struct;
       output.Alpha=alphaval;
       output.Omega=Omega;
       output.NumPerDriveCycle=numDrive;
       output.NumDriveCycles=numCycles;
       output.Final4CycleEnergy=mean(Edata(end-numDrive*4:end));
       output.Final4CycleDownStab=mean(omin_data(end-numDrive*4:end));
       output.Final4CycleUpStab=mean(omax_data(end-numDrive*4:end));
       output.Final4CycleStability=output.Final4CycleUpStab-output.Final4CycleUpStab;
       output.WF_data=psi_data(:,1:50:end);       
    end

    function theans=minOverlap(psi)
        rho=abs(psi).*abs(psi);
        theans=sum(rho.*rho_min);
    end

    function theans=maxOverlap(psi)
        rho=abs(psi).*abs(psi);
        theans=sum(rho.*rho_max); 
    end

    function E=staticEnergy(psi)
       E=abs(ctranspose(psi)*staticHmatrix*psi); 
    end
    
    function [val]=staticpotential(phi)
       val=-cos(phi)*a/(2*sqrt(a));
    end

    function [val]=relpotential(phi,theta)
       val=-cos(phi)*(a+b*cos(theta))/(2*sqrt(a)); 
    end

    function [psi_i] = makeGaussian(center)
    [psi_i]=arrayfun(@(phi) exp(-(sqrt(2*a)/8)*(phi-2*pi)^2),positionMesh);
    normal=sqrt(sum(abs(psi_i).*abs(psi_i)));
    
    [psi_i]=arrayfun(@(phi) exp(-(sqrt(2*a)/8)*(phi-center)^2),positionMesh);

    psi_i=psi_i/normal;
    psi_i=transpose(psi_i);
    end

    function [psi]=makeGaussianComb(state)         
       centers=-1:ceil(rangePhi/pi)+1;%Goes from -1 to N
       centers=2*pi*centers;%Now gos from -2pi to N*2pi
       centers=centers+pi*state;%if state=1 now goes from -pi to N*2pi-pi
        %if state=1 inverted
        %if state=0 nonverted
       psi=zeros(length(positionMesh),1);
       for index=1:length(centers);
           center=centers(index);
           thisState=makeGaussian(center);
           psi=psi+thisState;
       end
       
        normal=sqrt(sum(abs(psi).*abs(psi)));
        psi=psi/normal;        
        
    end
    
end

