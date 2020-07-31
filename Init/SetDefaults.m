function [Mat,Set ] = SetDefaults( Mat,Set)
% Sets default Values ofr structures Mat and Set

%%   Init\MainIncrements
if ~isfield(Set,'Profiler')
    Set.Profiler=false; % Matlab Profiler
end
if ~isfield(Set,'tol')
    Set.tol=1e-9; % Tolerance for Newton Raphson Convergence
end
if ~isfield(Set,'MaxIter')
    Set.MaxIter=20;
end
if ~isfield(Set,'BCcode')
    Set.BCcode=3; % =1: Incremental applied x-displacement at X=cont two boundaries. Stretching simulation.
                  % =-1: Same as 1, but only applied on the first step.
                  % =2: Applied force. Stretching simulation.
                  % =3: Fixed boundary. z for bottom, x and y for domain boundary. Wound healing simulation.
                  % =4: Fixed z for bottom, simulation of propulsion and friction (Extracellular Matrix) 
                  % =5: Free Boundary 
                  % =6: Fixed z for bottom with a region of free z bottom, 
                  %     The Region of wiht Free-Z defined by 
                  %     Set.ZFreeX=[X1 X2];


                  
                  
                  
                  
end 
if ~isfield(Set,'theta')
    Set.theta=0.5; % Value for time integration of viscous (resting length remodelling)
end
if ~isfield(Set,'tol')
    Set.tol=1e-9; % Tolerance for Newton Raphson Convergence
end
if ~isfield(Set,'MaxIter')
    Set.MaxIter=20;
end
% Maximum Number of StepHalvings
if ~isfield(Set,'StepHalvingMax')
    Set.StepHalvingMax=5;
end
% Maximum reduction of Set.dt during StepHalving. The minimum value of dt
% will be Set.dt/Set.StepHalvingMinDt,m with Set.dt the time-step when
% the first StepHalving started.
if ~isfield(Set,'StepHalvingMinDt')
    Set.StepHalvingMinDt=100; 
end
%% Ablation
if ~isfield(Set,'yRelaxation')
    Set.yRelaxation=true;
end 
if ~isfield(Set,'yRelaxationZ')
    Set.yRelaxationZ=false;
end

if ~isfield(Set,'AblationN')
    Set.AblationN=0;   % Number of ablated cells
    Set.EpsCTWE=0;   % Final value of applied contracitlity at vertices on top of wound edge
    Set.EpsCBWE=0;   % Final value of applied contracitlity at vertices on bottom of wound edge
    Set.EpsCLWE=0;   % Final value of applied contracitlity at vertices on lateral of wound edge
    Set.FixedBasal=0;
end 
if ~isfield(Set,'FixedBasal')
    Set.FixedBasal=0;
end 
if ~isfield(Set,'YmidWound')
    Set.YmidWound=1;         % =1 mid-plane vertices on the wound edge (In case Set.ModelTop=2 || =3)
                             % =0 without mid-plane vertices on the wound edge
end 
if ~isfield(Set,'AblationTimeStep')
    Set.AblationTimeStep=1;  % the time step at which ablation will take place. 
end 
if ~isfield(Set,'WRemodelThreshold') 
    Set.WRemodelThreshold=0.1;   % 0.1 Threshold for intercalation on the wound edge: it the retio 
end 
if ~isfield(Set,'RemodelEntangled') % Trey to remodel entangled elements
    Set.RemodelEntangled=true;
end
%%  Init\InitializeModel
if ~isfield(Set,'CellCentres')
    Set.CellCentres='';
end

if ~isfield(Set,'umPerPixel')
    Set.umPerPixel=1; % um per pixel of images. X coordinates in *.mat files are scales by this value
end
%% Geo
% Geo\MeshgenGrid
if ~isfield(Set,'Prisms')
    % =0; (default) Meshed with tets (Diagonalised nodal netwowrk),
    % =1 Meshed with triangular vertical prisms
    Set.Prisms=0;
end 
if ~isfield(Set,'nx') && ~isfield(Set,'CellCentres')
    warning('Set.nx not set, and no file CellCentres given. Set.nx=4')
    Set.nx=4;
end
if ~isfield(Set,'ny') && ~isfield(Set,'CellCentres')
    warning('Set.ny not set, and no file CellCentres given. Set.ny=4')
    Set.ny=4;
end
if ~isfield(Set,'nz')
    Set.nz=1;
end
if ~isfield(Set,'h')
    Set.h=1;
end
if ~isfield(Set,'XRand')
    Set.XRand=0.0; % Randomise initial positions. Recomended value ~0.2
end

%  Geo\MeshgenData
if ~isfield(Set,'ModelTop')
    Set.ModelTop=3; 
     % 1 = same top/bottom with no mid-plane vertices
     % 2 = same top/bottom with mid-plane vertices
     % 3 = different top/bottom with mid-plane vertices

end 

if ~isfield(Set,'RemodelExternal')
    Set.RemodelExternal=false; % false=Avoid remodelling on external boundary of domain
end
if ~isfield(Set,'RemodelDelta')
    Set.RemodelDelta=0.2; % Tolerance for graded Delaunay when Set.Remodel>0 
                          % RemodelDelta=0 standard Delaunay
                          % RemodelDelta>0 allows elongated cells.
                          % Recommended =0.2
end

if ~isfield(Set,'RemodelTolF')
    Set.RemodelTolF=5; % Tolerance for filtering boundary triangles in Delaunay Remodeling.
                         % r/R>TolF are filtered, r=cricumradius, R=inradius.
end

%%
if ~isfield(Set,'OutputVTK')
    Set.OutputVTK=1; % =1 write *.vtk files. =0 do not write files.
end


%% Material 
if ~isfield(Set,'lambdaV')
    Set.lambdaV=0; % Volume penalisation
end

if ~isfield(Mat,'Delay')
    Mat.Delay=0; % 0.4, no oscillation. 
               %  1= 2 oscil
               %  1.5: decaying oscillations. Vertex: milder oscillations
               %  2=sustained sudden oscillations, Vertex: instabilities.
               %  2.5= unstable at max 449 steps => tend=22.45)
               %  3=unstable latter oscil (max 332 stesps=>tend=16.6)

end 

if ~isfield(Mat,'D')
    Mat.D.k0=0;
end
if ~isfield(Mat,'V')
    Mat.V.k0=0;
end
if ~isfield(Mat.D,'k0')
    Mat.D.k0=0;
end
if ~isfield(Mat.D,'k')
    Mat.D.k=0;
end
if ~isfield(Mat.V,'k0')
    Mat.V.k0=0;
end
if ~isfield(Mat.V,'k')
    Mat.V.k=0;
end
if ~isfield(Mat.V,'gamma')
    Mat.V.gamma=0;
end

if ~isfield(Mat.D,'gamma')
    Mat.D.gamma=0;
end

if ~isfield(Mat.D,'EpsC')
    Mat.D.EpsC=0;% contractility Top
end
if ~isfield(Mat.V,'EpsCT')
    Mat.V.EpsCT=0;% contractility Top
end
if ~isfield(Mat.V,'EpsCB')
    Mat.V.EpsCB=0;% contractility Bottom
end
if ~isfield(Mat.V,'EpsCL')
    Mat.V.EpsCL=0;% contractility Lateral
end

%% Rheology
if ~isfield(Set,'Sparse')
    Set.Sparse=true; % Sparse Jacobian matrices 
end
if ~isfield(Set,'StrainBased')
    Set.StrainBased=1;% =1, Strain based: s=k*(l/L-1). =0, Displacement based: s=k*(l-L)
end
if ~isfield(Set,'Set.LRef0')
    Set.LRef0=true; % When true, and Set.StrainBased==false,  elastic branch uses sigma=k*l, instead of sigma=k*(l-L)
end

%% Substrate friction + Propulsion
%%% while using these options, make sure that the botttom vertices are not fully fixed 
if ~isfield(Set,'eta')
       Set.eta=0;        % friction viscosity
end 
if ~isfield(Set,'Propulsion')
    Set.Propulsion=false; 
%    if Set.Propulsion=true; it needs the following options 
%     % Region 1 with Propulsion forces Set.mu1  
%     Set.PropulsiveRegionX1=[0 round(Set.nx/5)];
%     Set.PropulsiveRegionY1=[0  Set.ny];
%     Set.mu1=[.1 0 0];       % velocity 
% 
%     % Region 2 with Propulsion forces Set.mu2  
%     Set.PropulsiveRegionX2=[round(Set.nx/1.33) Set.nx];
%     Set.PropulsiveRegionY2=[0  Set.ny];
%     Set.mu2=[-.1 0 0];       % velocity 
end 



end

