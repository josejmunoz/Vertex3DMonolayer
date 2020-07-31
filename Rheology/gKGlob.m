function [Cell,Ener,g,K,Set,L,Stress]=gKGlob(Ablated,C,Cv,Cell,Mat,L,L0,N,Stress,Set,i,T,X,Y,Yn,Ymid)

% Computes global residual and Jacobian

% Nodal network residual/Jacobian
[Ener.D.e,g,K,Stress.D.s]=gKNodal(Ablated,X,C,Mat.D,L.D.L,L0.D.L,Set);
% Vertex residual/Jacobian

[Ener.V.e,gy,Ky,Stress.V.s]=gKVertex(Ablated,Cell,Cv,Mat.V,L.V.L,L0.V.L,N,Y,Ymid,T,Set);





g=g+gy;
K=K+Ky;
% Volume residual/Jacobian
[gv,Kv,Cell]=gKVol(Cell,X,Y,Ablated.Yr,Set,T,N,Ymid);
g=g+gv;
K=K+Kv;


%---------------------- Propulsion Forces ------------------------------------
if Set.Propulsion
    [gp]=gPropulsion(Cell,Set);
     g=g-gp;
end

%---------------------- Viscous friction ------------------------------------

if Set.eta>0
    [gf,Kf]=gKFriction(Cell,N,T,Y,Yn,Set);
     K=K+Kf;
     g=g+gf;
end 