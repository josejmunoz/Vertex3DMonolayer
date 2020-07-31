function [Ld,ld,Mat]=DelayInit(C,Cv,Ln,Mat,X,Y,Set)
% Initiailsed delayed lengths
% INPUT:
% C  = connectivity of nodal lengths
% Cv = connectivity of vertex lengths
% X  = nodal positions
% Y  = vertex positions
% Set = structure with settings (see SetDefaults)
% OUTPUT:
% Mat = Material structure (see SetDefaults)
% Ld = set of previous resting lengths
% ld = set of previous apparent lengths
Mat.D.Delay=Mat.Delay;
Mat.V.Delay=Mat.Delay;
if Mat.Delay>0
    Mat.nDelay=ceil(Mat.Delay/Set.dt);
    Ld=cell(Mat.nDelay,1); % Delayed time-steps, resting lengths
    ld=cell(Mat.nDelay,1); % Delayed time-steps, current lengths
    Ld{1}=Ln;
    ld{1}=Lengths(C,Cv,X,Y);
    for i=2:Mat.nDelay
        Ld{i}=Ld{i-1};
        ld{i}=ld{i-1};
    end
else
    Mat.nDelay=0;
    Ld=cell(1);
    ld=cell(1);
    Ld{1}=Ln;
    ld{1}=Ln;
end
