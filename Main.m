% Main File for Vertex3DMonolayer
% See comments in InputData.m 
% Meaning of parameters are in Init/SetDefaults.m
%
% F Ioannou, J Munoz 20/11/2019
% Malik Dawi         20/07/2020
% Jose Munoz         25/07/2020 
%
% REFERENCES:
%
% Ioannou, F., Dawi, M.A., Tetley, R.J., Mao, Y., Munoz, J.J.
% Development of a new 3D hybrid model for epithelia morphogenesis
% Frontiers in Bioengineering and Biotechnology, Vol. 8, pp. 1-11, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
clc
tic
close
fclose('all');

% Input data with setting and Material properties
InputData % For running synthetic (cartesian) cell centres
%InputData99 % For reading file with 99 Cell Centres

%% Call to main function
if ~exist(strcat('.',Esc(),'Main.m'),'file')
    cd ..
    if ~exist(strcat('.',Esc(),'Main.m'),'file')
        error('Run Main from folder that contains file Main.m');
    end
end
addpath(strcat(pwd,Esc,'Init'));
[ C,Ener,X0,X,Y,Cell,Mat,Set,Stress,Ablated] = MainIncrements(Mat, Set );