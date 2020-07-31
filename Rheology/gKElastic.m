%% Elemental residual for elastic element
function [ge,ke,s,ener]=gKElastic(k,Mat,epsC,L,X,Set,wound,branch)
% INPUT:
%  k  = elemental stiffness (may be idfferent from the general one in Mat)
%  X  = Current Nodal/vertex coordinates
%  X0 = Initial Nodal/vertex coordinates
%  branch =1. Purely elastic. Stress=k(l(L-1)
%         =2. Contractile viscoelastic. Stress=k(l/L-1) and dot L=k*(l-L(1+epsC))
dim=size(X,2);
x1=X(1,:);
x2=X(2,:);
l=norm(x2-x1);
if Set.StrainBased
    if abs(L)<eps
        strain=0;
    else
        strain=l/L-1;
    end
else
    if Set.LRef0
        strain=l;
    else
        strain=l-L;
    end
end
ge=zeros(2*dim,1);

if branch==1
    epsCW=0;
%             gamma=Mat.gamma;
%         [epsCW,~]=CheckAblated(wound,Set,gamma);

    dLdl=0;
else
    gamma=Mat.gamma;
    [epsCW,gamma]=CheckAblated(wound,Set,gamma);
    if abs(Mat.Delay)>eps
        dLdl=0;
    else
        dLdl=Set.dt*gamma*Set.theta/(1+Set.dt*gamma*Set.theta);  %%% %%% %%%
    end
end
epsC=epsC+epsCW; % Final contracility=Material(undounded) + DueToWounding
% if branch==2   %%%
%     epsC=0;    %%%
% end            %%%
strainT=strain+epsC; % Add contracility
if abs(l)<eps
    ge(1:dim)=0;
    ge(dim+1:2*dim)=0;
else
    ge(1:dim)=k*strainT/l*(x1-x2); % D_x Phi^e
    ge(dim+1:2*dim)=k*strainT/l*(x2-x1); % =-ge(1:dim)
end
s=k*strainT;
ener=0.5*k*[strain^2 epsC^2 epsCW^2];
%%%%
n=size(X,1);
ke=zeros(3*n);
e=zeros(dim,2);
e(:,1)=(x2-x1)/l;
e(:,2)=(x1-x2)/l;
for i=1:n
    idof=(i-1)*dim+1:i*dim;
    for j=1:n
        jdof=(j-1)*dim+1:j*dim;
        if abs(l)<eps || (Set.StrainBased && abs(L)<eps)
            ke(idof,jdof)=(-1)^(i+j)*k*eye(dim);
        else
            if Set.StrainBased
                ke(idof,jdof)=(-1)^(i+j)*(k*strainT/l)*eye(dim)+... % Gamma
                    k*(1/L-strainT/l-l/L^2*dLdl)*e(:,i)*e(:,j)';
            else
                ke(idof,jdof)=(-1)^(i+j)*(k*strainT/l)*eye(dim)+... % Gamma
                    k*(1-strainT/l-dLdl)*e(:,i)*e(:,j)';
            end
        end
    end
end
end

