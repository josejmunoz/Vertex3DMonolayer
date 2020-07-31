function [epsC,gamma]=CheckAblated(wound,Set,gamma)
% Updated to be compatible with Rheology\Ablation.m 
% INPUT
% wound= 0 healthy element 
%      =1: Ablated, nodal top
%      =2: Ablated, nodal bottom
%      =3: Ablated, nodal vertical
%      = 4 ablated bottom element
%      = 5 ablated top element
%      = 6 ablated lateral (all type of Lateral)
%      = 7 wounded bottom element  (on the wound edge)
%      = 8 wounded top element     (on the wound edge)
%      = 9 wounded Lateral element (on the wound edge)
% OUTPUT
% epsC :   = value of contractility on wounded elements
% Contractility
epsC=0;
switch(wound)
    case(1)
        epsC=Set.EpsCTW_N*Set.t(Set.iIncr)/Set.t(Set.Nincr);
    case(2)
        epsC=Set.EpsCBW_N*Set.t(Set.iIncr)/Set.t(Set.Nincr);
    case(3)
        epsC=Set.EpsCLW_N*Set.t(Set.iIncr)/Set.t(Set.Nincr);
    case(4)
        epsC=Set.EpsCBW*Set.t(Set.iIncr)/Set.t(Set.Nincr);
    case(5)
        epsC=Set.EpsCTW*Set.t(Set.iIncr)/Set.t(Set.Nincr);
    case(6)
        epsC=Set.EpsCLW*Set.t(Set.iIncr)/Set.t(Set.Nincr);
    case(7)
%         epsC=Set.EpsCBWE*Set.iIncr/Set.Nincr;
          epsC=Set.EpsCWBot(Set.iIncr); 
    case(8)
%         epsC=Set.EpsCTWE*Set.iIncr/Set.Nincr;
          epsC=Set.EpsCWTop(Set.iIncr); 

    case(9)
%         epsC=Set.EpsCLWE*Set.iIncr/Set.Nincr;
          epsC=Set.EpsCWLat(Set.iIncr); 

end
% Vsicosity: decrease visocisty (increase gamma) if inside wound
if nargin>2 && wound<7 && wound>0
    gamma=gamma*Set.gammaFact;
end
end