%%   Init\MainIncrements
Set.BCcode=3; % =1: Incremental applied x-displacement at X=cont two boundaries. Stretching simulation.
              % =-1: Same as 1, but only applied on the first step.
              % =2: Applied force. Stretching simulation.
              % =3: Fixed boundary nodes. z for bottom, x and y for domain boundary. Wound healing simulation.    
Set.yRelaxation=true;
Set.tend=150; % 40s
Set.dt=1; % Time-step size (in wound healing applied during closure)
Set.t=[0:0.6:6 6+Set.dt:Set.dt:Set.tend]; % Time history. If not declared, Set.t=0:Set.dt:Set.tend
% Set.t=[0:0.6:6 12:3:72]; % Filippos

Set.CellCentres='CellCentres99.mat';%'Xtb.mat';% CellCentres.mat=62 cells. CellCentres1.mat=153 cells, CellCentres221.mat=221 cells, ...

Set.nz=1;
Set.h=20; % Height in um
Set.umPerPixel=0.0527; % Scaling of X coordinates. Units in *.mat files with cell centres ar in pixel units.

Set.ModelTop=1; 
     % 1 = same top/bottom with no mid-plane vertices
     % 2 = same top/bottom with mid-plane vertices
     % 3 = different top/bottom with mid-plane vertices

Set.RemodelDelta=5;%0.05; % Tolerance for graded Delaunay when Set.Remodel>0 
                      % RemodelDelta=0 standard Delaunay
                      % RemodelDelta>0 allows elongated cells.
                      % Recommended =0.2


Set.RemodelTolF=5; % Tolerance for filtering boundary triangles in Delaunay Remodeling.
                     % r/R>TolF are filtered, r=cricumradius, R=inradius.

%% Material 
% volume 
Set.lambdaV=20; % Volume penalisation

% Delaunay
Mat.D.k0=0.3; % Stiffness Delaunay Spring branch
Mat.D.k=0.0; % Stiffness Delaunay Spring branch
Mat.D.EpsC=0.3; % contractility (no distinction top, bottom ,lateral on nodal)
% Vertex
Mat.V.k0=0.05; % Stiffness Vertices spring branch
Mat.V.k=1.0; % 1; Stiffness Vertices spring branch
Mat.V.gamma=0.2; % 0.8
Mat.V.EpsCT=1.3; %0.2 contractility Top
Mat.V.EpsCB=1.3; %0.2 contractility Bottom     
Mat.V.EpsCL=.01; % 0.05 contractility Laterals
%% Ablation
Set.AblationN=5;   % Number of ablated cells
Set.EcTypeBot=2;       % Type of The contractility  profile. See 'help ContractilityInit'
                       % =1 Step function (needed parameter Set.StartTimeEcBot, Set.EpsCBWE)
                       % =2 Hat  function  (needed parameter Set.StartTimeEcBot, Set.PeakTimeEcBot, Set.EndTimeEcBot, Set.EpsCBWE)
Set.EpsCBWE=0.5;    % value of applied contracitlity at vertices on bottom of wound edge;
Set.StartTimeEcBot=9+Set.dt;  % The time at which the contractility start to be applied
Set.PeakTimeEcBot=16;   % The time at which the contractility reaches the prescribed value (hat function)
Set.EndTimeEcBot=700;    % The time at whihc the contractility reduced to zero (hat  function
   
Set.EcTypeTop=2;  
Set.EpsCTWE=2.3;
Set.StartTimeEcTop=9+Set.dt;  
Set.PeakTimeEcTop=16; % 6+6+Set.dt;   
Set.EndTimeEcTop=500;    
    
Set.EcTypeLat=1;  
Set.EpsCLWE=0.01;    
Set.StartTimeEcLat=13;  
Set.PeakTimeEcLat=13;   
Set.EndTimeEcLat=0; 
    
Set.FixedBasal=1;
Set.AblationTimeStep=1;      % the time step at which ablation will take place. 
Set.YmidWound=1;             % =1 mid-plane vertices on the wound edge
                             % =0 without mid-plane vertices on the wound edge
Set.WRemodelThreshold=0.15;   % 0.1 Threshold for intercalation on the wound edge: it is the maximum allowed aspect ratio

Set.OutputVTK=true;
Set.MaxIter=15;
Set.StepHalvingMax=3;
Set.StepHalvingMinDt=10; % Maximum redution of dt
Set.RemodelEntangled=false;

