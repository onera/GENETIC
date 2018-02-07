% GENETIC.BENCH.LOAD enables to load benchmark.
%
% Syntax
%  b = genetic.bench.load(name)
%  b = genetic.bench.load(name, xDim)
%
% Parameters
%  name : name of the benchmark
%  xDim : dimension of the search space (when it is not given by the
%  bench itself)
%
% Example
%  xDim                 = 10
%  ben                  = genetic.bench('ackley', xDim);
%  opt.verbosity        = 1;
%  [xopt, fopt, info]   = genetic.min({ben.f, ben.fDim}, ben.xDim, 'cmaes', opt)
%
%  Note that the interface genetic.min enables to call directly a benchmark:
%
%  opt.verbosity        = 1;
%  [xopt, fopt, info]   = genetic.min('ackley', 10, 'cmaes', opt);
%

% -------------------------------------------------------------------------

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
function bench = load(benchName, n ,quiet)
bench = [];

if nargin < 3
   quiet = false;
end
if nargin >= 2 && ~isempty(n)
   xDim = n;
else
   xDim = [];
end
if nargin >= 1
   [benchExists, path]  = genetic.bench.checkAvailability(benchName);
   if benchExists
      bench             = loadData(benchName, path, xDim, quiet);
   else
      error('The benchmark ''%s'' does not appear to exist. Please verify the spelling.\n',self.name);
   end
end
end

function bench = loadData(benchName, path, xDim, quiet)
bench = [];       
if isempty(strfind(path, ['+bench',filesep,'+data',filesep,'+usr']))
   % genetic benchs
   pattern     = ['+genetic',filesep,'+bench',filesep,'+data'];
   tmpBench    = path(strfind(path,pattern)+length(pattern)+1:end);
   tmpBench    = strrep(tmpBench,filesep(),'.');
   tmpBench    = strrep(tmpBench,'+','');
   benchCall   = ['genetic.bench.data.',tmpBench];
   data        = feval(benchCall);
else
   % user defined bench
   data        = importdata([path,'.mat']);
end
if isempty(data)
   if isempty(xDim)
      if ~quiet
         fprintf('The dimension of the benchmark ''%s'' must be specified.\n',benchName);
      end
      return
   else
      data        = feval(benchCall, xDim);
      data.xDim   = xDim;
   end
else
   if ~isempty(xDim)  && xDim ~= data.xDim
      if ~quiet
         fprintf('The dimension of the benchmark ''%s'' is fixed to %d. Discarding the provided value of %d\n',benchName, data.xDim, xDim);
      end
   end
end
bench       = genetic.bench.getEmpty();
data.name   = benchName;
bench       = genetic.tools.fillStructure(bench, data);
end