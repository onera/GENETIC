% Copyright 2018 ONERA
%
% This file is part of the GENETIC project.
%
% GENETIC is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License version 3 as
% published by the Free Software Foundation.
%
% GENETIC is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with GENETIC.  If not, see <https://www.gnu.org/licenses/lgpl-3.0>.
%
function [fun, n, cst] = toCall(bench)
% Creates the elements required for calling genetic from a benchmark
% structure
fun   = {bench.f, bench.fDim};
n     = bench.xDim;
cst   = struct('A',bench.A,'b',bench.b,'Aeq',bench.Aeq,'beq',bench.beq,'ceq',bench.ceq,'c',bench.c,'xMin',bench.bounds(:,1),'xMax', bench.bounds(:,2));
end