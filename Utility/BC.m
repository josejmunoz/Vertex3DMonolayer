function [gext,dofP,U0 ] = BC(Set,X,xExternal)
%BC set boundary conditions
%   BCcode     =1: Applied displacement. Stretching simulation.
%              =2: Applied force. Stretching simulation.
%              =3: Fixed boundary. z for bottom, x and y for domain boundary. Wound healing simulation.
dim=size(X,2);
nodes=Set.nodes;
nvert=Set.nvert;
BCcode=Set.BCcode;
if Set.yRelaxation
    gext=zeros((nodes+nvert)*dim,1);
else
    gext=zeros(nodes*dim,1);
end
%
if abs(BCcode)==1  % Applied x displacements
    u0=Set.u0; 
    xmin=min(X(:,1));
    xmax=max(X(:,1));
    nodLeft=find(X(:,1)==xmin);
    nodRight=find(X(:,1)==xmax);
    nLeft=size(nodLeft,1);
    nRight=size(nodRight,1);
    U0=[nodLeft ones(nLeft,1) zeros(nLeft,1)
        nodLeft 2*ones(nLeft,1) zeros(nLeft,1)
        nodLeft 3*ones(nLeft,1) zeros(nLeft,1)
        nodRight ones(nRight,1) u0*ones(nRight,1)
        nodRight 2*ones(nRight,1) zeros(nRight,1)
        nodRight 3*ones(nRight,1) zeros(nRight,1)        ];
elseif BCcode==2 % Applied force
    if ~isfield(Set,'nx') && ~isfield(Set,'ny')
        error('when using BCcode=2 (applied force), values for Set.nx and Set.ny are needed.')
    end
    nx=Set.nx;
    ny=Set.ny;
    U0=zeros((ny+1)*dim,3);
    F=zeros((ny+1)*dim,3);
    t=0;
    tF=0;
    for i=1:ny+1
        ii=i;
        for j=1:dim % Left node
            t=t+1;
            U0(t,:)=[ii,j,0];
        end
        ii=i+nx*(ny+1);
        tF=tF+1;
        F(tF,:)=[ii,1,1];
        for j=2:dim % Left node
            tF=tF+1;
            F(tF,:)=[ii,j,0];
        end
    end
    gext((F(:,1)-1)*dim+F(:,2))=F(:,3);
    dofP=(U0(:,1)-1)*dim+U0(:,2); % Prescribed dof
elseif BCcode==3 % Fixed boundary  
    zmin=min(X(:,3));
    nodBot=find(X(:,3)==zmin);
    nBot=size(nodBot,1);
    if isempty(xExternal)
        xmin=min(X(:,1));
        xmax=max(X(:,1));
        ymin=min(X(:,2));
        ymax=max(X(:,2));
        nodLeft=find(X(:,1)==xmin);
        nodRight=find(X(:,1)==xmax);
        nodDown=find(X(:,2)==ymin);
        nodUp=find(X(:,2)==ymax);
        nLeft=size(nodLeft,1);
        nRight=size(nodRight,1);
        nUp=size(nodUp,1);
        nDown=size(nodDown,1);
        U0=[nodLeft   ones(nLeft,1)  zeros(nLeft,1)
            nodRight  ones(nRight,1) zeros(nRight,1)
            nodUp    2*ones(nUp,1)   zeros(nUp,1)
            nodDown  2*ones(nDown,1) zeros(nDown,1)
            nodBot   3*ones(nBot,1)  zeros(nBot,1)  ];
    else
        nodes=size(X,1)/2;
        xExternal=[xExternal xExternal+nodes];
        nExt=length(xExternal);
        U0=[xExternal' ones(nExt,1) zeros(nExt,1)
            xExternal' 2*ones(nExt,1) zeros(nExt,1)
            nodBot   3*ones(nBot,1)  zeros(nBot,1)  ];
    end
    dofP=(U0(:,1)-1)*dim+U0(:,2); % Prescribed dof
elseif BCcode==4 % Free Boundary (Fixed z for bottom)
    zmin=min(X(:,3));
    nodBot=find(X(:,3)==zmin);
    nBot=size(nodBot,1);
    U0=[nodBot   3*ones(nBot,1)  zeros(nBot,1)  ];
    dofP=(U0(:,1)-1)*dim+U0(:,2); % Prescribed dof
%     U0=[nan nan nan];
%     dofP=[]; % Prescribed dof
elseif BCcode==5 % Free Boundary 
    U0=[nan nan nan];
    dofP=[]; % Prescribed dof
elseif BCcode==6 % Free Boundary (Fixed z for bottom but free in some region defined by (Set.ZFreeX) 
    ID=1:size(X,1);
    nodBot=ID(X(:,1)<Set.ZFreeX(1) | X(:,1)>Set.ZFreeX(2))';
    nBot=size(nodBot,1);
    U0=[nodBot   2*ones(nBot,1)  zeros(nBot,1)  
        nodBot   3*ones(nBot,1)  zeros(nBot,1)];
    dofP=(U0(:,1)-1)*dim+U0(:,2); % Prescribed dof

else
    error('Wrong BCcode')
end
end

function U0=U0Manual()
%%% X-AXIS DISPLACEMENT
U0=[1,1,0
    2,1,0
    3,1,0
    4,1,0
    5,1,0
    6,1,0
    7,1,0
    8,1,0
    9,1,0
    10,1,0
    11,1,0
    12,1,0
    13,1,2  %
    14,1,2  %
    15,1,2  %
    16,1,2  %
    17,1,0
    18,1,0
    19,1,0
    20,1,0
    21,1,0
    22,1,0
    23,1,0
    24,1,0
    25,1,0
    26,1,0
    27,1,0
    28,1,0
    29,1,2   %
    30,1,2   %
    31,1,2   %
    32,1,2   %
    1,2,0
    2,2,0
    3,2,0
    4,2,0
    5,2,0
    6,2,0
    7,2,0
    8,2,0
    9,2,0
    10,2,0
    11,2,0
    12,2,0
    13,2,0  %
    14,2,0  %
    15,2,0  %
    16,2,0  %
    17,2,0
    18,2,0
    19,2,0
    20,2,0
    21,2,0
    22,2,0
    23,2,0
    24,2,0
    25,2,0
    26,2,0
    27,2,0
    28,2,0
    29,2,0  %
    30,2,0  %
    31,2,0   %
    32,2,0   %
    1,3,0
    2,3,0
    3,3,0
    4,3,0
    5,3,0
    6,3,0
    7,3,0
    8,3,0
    9,3,0
    10,3,0
    11,3,0
    12,3,0
    13,3,0
    14,3,0
    15,3,0
    16,3,0
    17,3,0
    18,3,0
    19,3,0
    20,3,0
    21,3,0
    22,3,0
    23,3,0
    24,3,0
    25,3,0
    26,3,0
    27,3,0
    28,3,0
    29,3,0
    30,3,0
    31,3,0
    32,3,0];
end