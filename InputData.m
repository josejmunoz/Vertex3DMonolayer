% Inputa Data for 3D analysis of  wound healing on monolayers
%
% v4.9
% 20/7/2020
% REFERENCE:
% 
% Ioannou, F., Dawi, M.A., Tetley, R.J., Mao, Y., Munoz, J.J.
% Development of a new 3D hybrid model for epithelia morphogenesis
% Frontiers in Bioengineering and Biotechnology, Vol. 8, pp. 1-11, 2020
% 

Set.yRelaxation=true; % trueBottom verties no not follow 
Set.tend=15;          % Final time
Set.dt=1;             % Time-step size (in wound healing applied during closure)
% May be declared as a vector:
% Set.t=[0:0.6:6 6+Set.dt:Set.dt:Set.tend]; % Time history. If not declared, Set.t=0:Set.dt:Set.tend
% Set.t=[0:0.6:6 12:3:72]; % Filippos

% GEOMETRY:
% For experimental location of Cell Centres:
% Set.CellCentres='CellCentres106.mat';%'Xtb.mat';% CellCentres.mat=62 cells. CellCentres1.mat=153 cells, CellCentres221.mat=221 cells, ...

% Geometrical propoerties if not experimental lovation
Set.nx=8;
Set.ny=8;
Set.nz=1;
Set.h=3; % Height in um
Set.umPerPixel=0.0527; % Scaling of X coordinates. Units in *.mat files with cell centres ar in pixel units.

% BOUNDARY CONDITIONS:
Set.BCcode=3; % =1: Incremental applied x-displacement at X=cont two boundaries. Stretching simulation.
              % =-1: Same as 1, but only applied on the first step.
              % =2: Applied force. Stretching simulation.
              % =3: Fixed boundary. z for bottom, x and y for domain boundary. Wound healing simulation.
              % =4: Fixed z for bottom, simulation of propulsion and friction (Extracellular Matrix) 
              % =5: Free Boundary 
              % =6: Fixed z for bottom with a region of free z bottom, 
              %     The Region of wiht Free-Z defined by 
              %     Set.ZFreeX=[X1 X2];    
                Set.ZFreeX=[round(Set.nx/5) round(Set.nx/1.2500)];
Set.ModelTop=1; 
     % 1 = same top/bottom with no mid-plane vertices
     % 2 = same top/bottom with mid-plane vertices
     % 3 = different top/bottom with mid-plane vertices

Set.RemodelDelta=.2;%0.05; % Tolerance for graded Delaunay when Set.Remodel>0 
                      % RemodelDelta=0 standard Delaunay
                      % RemodelDelta>0 allows elongated cells.
                      % Recommended =0.2


Set.RemodelTolF=50; % Tolerance for filtering boundary triangles in Delaunay Remodeling.
                     % r/R>TolF are filtered, r=cricumradius, R=inradius.

%% MATERIAL PROPERTIES 
% volume 
Set.lambdaV=20; % Volume penalisation

% Delaunay
Mat.D.k0=0.5; % Stiffness Delaunay Spring branch
Mat.D.k=0.0; % Stiffness Delaunay Spring branch
Mat.D.EpsC=0.0; % contractility (no distinction top, bottom ,lateral on nodal)
% Vertex
Mat.V.k0=0.5; % Stiffness Vertices spring branch
Mat.V.k=1.0; % 1; Stiffness Vertices spring branch
Mat.V.gamma=0.2; % 0.8
Mat.V.EpsCT=0.10; %0.2 contractility Top
Mat.V.EpsCB=0.10; %0.2 contractility Bottom     
Mat.V.EpsCL=0.40; % 0.05 contractility Laterals

%% Ablation
Set.AblationN=1;   % Number of ablated cells
Set.EcTypeBot=2;       % Type of The contractility  profile. See 'help ContractilityInit'
                       % =1 Step function (needed parameter Set.StartTimeEcBot, Set.EpsCBWE)
                       % =2 Hat  function  (needed parameter Set.StartTimeEcBot, Set.PeakTimeEcBot, Set.EndTimeEcBot, Set.EpsCBWE)
Set.EpsCBWE=0;    % value of applied contracitlity at vertices on bottom of wound edge;
Set.StartTimeEcBot=6;  % The time at which the contractility start to be applied
Set.PeakTimeEcBot=18;   % The time at which the contractility reaches the prescribed value (hat function)
Set.EndTimeEcBot=24;    % The time at whihc the contractility reduced to zero (hat  function
   
Set.EcTypeTop=2;  
Set.EpsCTWE=2.2;
Set.StartTimeEcTop=1;  
Set.PeakTimeEcTop=16; % 6+6+Set.dt;   
Set.EndTimeEcTop=800;    
    
Set.EcTypeLat=1;  
Set.EpsCLWE=0.01;    
Set.StartTimeEcLat=13;  
Set.PeakTimeEcLat=13;   
Set.EndTimeEcLat=0; 
    
Set.FixedBasal=1;
Set.AblationTimeStep=1;      % the time step at which ablation will take place. 
Set.YmidWound=1;             % =1 mid-plane vertices on the wound edge
                             % =0 without mid-plane vertices on the wound edge
Set.WRemodelThreshold=0.1;   % 0.1 Threshold for intercalation on the wound edge: it is the maximum allowed aspect ratio

Set.OutputVTK=true;
Set.MaxIter=15;
Set.StepHalvingMax=3;

%% Substrate friction + Propulsion
%%% while using these options, make sure that the botttom vertices are not fully fixed 

Set.eta=.0;        % friction viscosity

Set.Propulsion=false; 
% Region 1 with Propulsion forces Set.mu1  
Set.PropulsiveRegionX1=[0 round(Set.nx/5)];
Set.PropulsiveRegionY1=[0  Set.ny];
Set.mu1=[.2 0 0];       % velocity 

% Region 2 with Propulsion forces Set.mu2  
Set.PropulsiveRegionX2=[round(Set.nx/1.33) Set.nx];
Set.PropulsiveRegionY2=[0  Set.ny];
Set.mu2=[-.2 0 0];       % velocity 