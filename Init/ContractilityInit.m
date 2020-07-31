function [Set]=ContractilityInit(Set)
%% Initialises contractiity profile for Bottom/Top/Lateral bars 
% INPUT:
% Set.EcTypeXX = 1 Use Step function for XX=Top,Bot,Lat
%              = 2 Use Hat function for XX=Top,Bot,Lat
% StartTimeEcXX   
% PeakTimeEcXX     Timing settings (see diagram below)
% EndTimeEcXX
% Eps
%  STEP FUNCTION      HAT FUNCTION
%    ^    EpsC           ^           EpsC       
%    |   +-------        |           / |\         
%    |   |               |         /   | \        
%    |   |               |       /     |  \       
%    |   |               |     /       |   \      
%    +---+-------->      +---+---------+----+---> 
%       t0                   t0       t1   t2     
%   t0=StartTimeEc
%   t1=PeakTimeEc
%   t2=EndTimeEc
% OUTPUT:
% Set.EpsCWXX
%%
if     Set.EcTypeBot==1 %------------------ Step function
         EpsCWBot=zeros(Set.Nincr,1);
         nt1=ceil(Set.StartTimeEcBot/Set.dt);
         EpsCWBot(nt1:end)=Set.EpsCBWE;
elseif Set.EcTypeBot==2 %------------------ Hat function 
         [EpsCWBot]=HatProfile(Set.StartTimeEcBot,....
                              Set.PeakTimeEcBot,...
                              Set.EndTimeEcBot,...
                              Set.EpsCBWE,Set.t);
end

if     Set.EcTypeTop==1 %------------------ Step function
         EpsCWTop=zeros(Set.Nincr,1);
         nt1=ceil(Set.StartTimeEcTop/Set.dt);
         EpsCWTop(nt1:end)=Set.EpsCBWE;
elseif Set.EcTypeTop==2 %------------------ Hat function 
         [EpsCWTop]=HatProfile(Set.StartTimeEcTop,....
                              Set.PeakTimeEcTop,...
                              Set.EndTimeEcTop,...
                              Set.EpsCTWE,Set.t);
end 
% figure(6)
% plot(1:Set.tend,EpsCWTop)
if     Set.EcTypeLat==1 %------------------ Step function
         EpsCWLat=zeros(Set.Nincr,1);
         nt1=ceil(Set.StartTimeEcLat/Set.dt);
         EpsCWLat(nt1:end)= Set.EpsCLWE;
elseif Set.EcTypeLat==2 %------------------ Hat function 
         [EpsCWLat]=HatProfile(Set.StartTimeEcLat,....
                              Set.PeakTimeEcLat,...
                              Set.EndTimeEcLat,...
                              Set.EpsCLWE,Set.t);
end 
% figure(7)
% plot(1:Set.tend,EpsCWLat)
Set.EpsCWTop=EpsCWTop;
Set.EpsCWBot=EpsCWBot;
Set.EpsCWLat=EpsCWLat;
end

%% Local functions
function EpsC=HatProfile(StartTime,PeakTime,EndTime,epsC,t) 
% Returns array with Hat Porfile of contactility
EpsC=zeros(length(t),1);
t0=StartTime;
t1=PeakTime;
t2=EndTime;
nt1=find(t>=StartTime,1);
if isempty(nt1)
    nt1=length(t);
end
%         p1=1:nt1;
nt2=find(t>=PeakTime,1);
if isempty(nt2)
    nt2=length(t);
end
p2=nt1:nt2;
nt3=find(t>=EndTime,1);
if isempty(nt3)
    nt3=length(t);
end
p3=nt2:nt3;
if t1>t0
    slop1=epsC/(t1-t0);
    EpsC(p2)=(t(p2)-t0).*slop1;
end
if t2>t1
    slop2=epsC/(t2-t1);
    EpsC(p3)=epsC-(t(p3)-t1).*slop2;
end
end
% Set.AblationN=1;   % Number of ablated cells
% Set.EpsCTWE=0;   % Final value of applied contracitlity at vertices on top of wound edge
% Set.EpsCBWE=0;               % Final value of applied contracitlity at vertices on bottom of wound edge
% Set.EpsCLWE=1;               % Final value of applied contracitlity at vertices on lateral of wound edge
% Set.FixedBasal=1;
% Set.AblationTimeStep=3;      % the time step at which ablation will take place. 
% Set.YmidWound=1;             % =1 mid-plane vertices on the wound edge
%                              % =0 without mid-plane vertices on the wound edge
% Set.WRemodelThreshold=0.1;   % 0.1 Threshold for intercalation on the wound edge: it the retio 