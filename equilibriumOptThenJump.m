function [height,state,o]=equilibriumOptThenJump(fi,tstart,P,varargin)
% function [height,state,o]=equilibriumOptThenJump(tstart,fi,P,varargin)
% This function computes the jump height. squat jump, no counter-movements.
% INPUTS:
% tstart = onset timing of the 6 muscles.
% fi = initial angles of the toe, ankle, knee and hip joints.
% P = structure of model parameters. this is currently done by loading a
% table of data from getJumperParams (returns a matrix) and then
% getJumperStruct(inMatrix)


%% optimization chunk.
fprintf('Computing initial starting gamma, stim, and muscle lengths for equilibrium.\n');
fi=fi(:);
[clcerel,cgamma,istim,ctor,~]=equilibriumOpt(fi',P);

istate = [fi(:)',zeros(1,4),zeros(1,2),zeros(1,2),clcerel(:)',cgamma(:)'];
if isempty(varargin)
    [height,state,o]=jumpMuscle(istate,istim,tstart,P);
else
    [height,state,o]=jumpMuscle(istate,istim,tstart,P,varargin{:});
end;
