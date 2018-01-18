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
function create(stru)
% GENETIC.BENCH.CREATE creates a benchmark from the data in the
% structure. The benchmark is stored in the folder +usr of +data
% and can be loaded as the other benchmarks.
% The dimension of the space cannot be adjusted though.
%
% Example
%  ben.f    = @(x) norm(x);
%  ben.xDim = 10;
%  ben.fopt = 0;
%  ben.xopt = zeros(10,1);
%  ben.name = 'norm10';
%  genetic.bench.create(ben);
%  clear all
%  ben      = genetic.bench.load('norm10')
if genetic.bench.checkAvailability(stru.name)
   fprintf('The benchmark ''%s'' already exists. Cannot create it.\n',stru.name);
   return
end
info        = what('+genetic');
geneticPath = info.path;
target      = [geneticPath, filesep, '+bench',filesep,'+data', filesep, '+usr', filesep, stru.name,'.mat'];
%
emptyBench  = genetic.bench.getEmpty();
bench       = genetic.tools.fillStructure(emptyBench, stru);
%
save(target,'bench');
end
