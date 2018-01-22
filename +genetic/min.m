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
function [xopt, fopt, info] = min(fun, x, method, varargin)
% GENETIC.MIN is the interface for minimisation in Genetic.
%
% Syntax
%  [xopt, fopt, info] = GENETIC.MIN(fun, n, method)
%  [xopt, fopt, info] = GENETIC.MIN(fun, n, method, options)
%  [xopt, fopt, info] = GENETIC.MIN(fun, n, method, constraints, options)
%  [xopt, fopt, info] = GENETIC.MIN({fun, nobj}, n, method, constraints, options)
%
% Parameters
%  fun         : objective function to minimise (handle function)
%  nobj        : number of objective (integer in a cell with the function)
%  n           : dimension of the optimisation problem or initial point (integer
%                or array)
%  method      : name of the method to use (string)
%  options     : parameters for the optimisation process (structure)
%  constraints : constraints of the optimisation problem (structure)
%  xopt        : best solution or front of solutions (array)
%  fopt        : best objective function found or Pareto front (array)
%  info        : informations about the optimisation process (structure)
%
% Description
%
% Example
%  % Problem formulation
%  n = 3;
%  A = rand(n,n);
%  b = rand(n,1);
%  f = @(x) norm(A*x - b,inf);
%  % Optimisation
%  options.maxFunEval   = 5000;
%  options.verbosity    = 1;
%  [xopt, fopt, info]   = genetic.min(f, n, 'cmaes', options);
%
% See also
%  <a href="matlab: help genetic.methods">methods</a>
%  <a href="matlab: help genetic.parameters">parameters</a>
%

%
tStart      = tic();
constraints = [];
options     = [];
%%
% =========================================================================
% Input arguments handling
% =========================================================================
%
% Initial point / dimension of the problem ................................
%
if ~isempty(x)
   if length(x) > 1 || rem(x,1) ~= 0
      % if x is a vector or a floating point value, then it is considered as
      % the initial point and the dimension of the problem its row size.
      x0    = x;
      xDim  = size(x,1);
   else
      % otherwise, there is no initial point and the dimension of the problem
      % is directly x
      x0    = [];
      xDim  = x;
   end
else
   % The dimension/initial point can be empty for some benchmarks (as it is
   % fixed by the benchmark itself)
   if ~ischar(fun)
      error('Problem dimension/initial point cannot be empty')
   end
   xDim  = [];
   x0    = [];
end
%
% Objective function/ multi-objective case ................................
%
if ischar(fun)
   % Direct call to a benchmark from genetic
   bench       = genetic.bench.load(fun,xDim);
   if isempty(bench.xDim)
      error('The dimension of the benchmark must be specified')
   end
   fun         = {bench.f, bench.fDim};
   constraints = struct('A'   , bench.A,...
                        'b'   , bench.b,...
                        'Aeq' , bench.Aeq,...
                        'beq' , bench.beq,...
                        'c'   , bench.c,...
                        'ceq' , bench.ceq,...
                        'xMin', bench.bounds(:,1),...
                        'xMax', bench.bounds(:,2));
%
   xDim        = bench.xDim;
end
if iscell(fun)
   % if fun is a cell, then it is assumed that the first element is the
   % objective function (handle) while the second element is the dimension
   % of the functions
   nobj  = fun{2};
   fun   = fun{1};
else
   % otherwise, the problem is assumed to be mono-objective
   nobj  = 1;
end
%
simulator            = genetic.simulator(fun, xDim, nobj);
%
% Method to be used .......................................................
%
if nobj == 1
   % Mono-objective methods
   availableMethods = genetic.tools.listMethods('mo');
else
   % Multi-objective methods
   availableMethods = genetic.tools.listMethods('mu');
end
%
method = lower(method);
if isempty(intersect(method, availableMethods))
   error('Method ''%s'' is not available',method);
end
%
% Constraints and options .................................................
%
[constraints, options] = genetic.tools.separateConstraintsAndOptions(varargin);
%
constraints          = genetic.constraints(constraints);
% Assignin verbosity to genetic
genetic.tools.params('verbosity', 0);
if isfield(options, 'verbosity')
   genetic.tools.params('verbosity', options.verbosity);
   options = rmfield(options, 'verbosity');
end
% Some simulator options
simOptions = {'gradObj'};
for i = 1:length(simOptions)
   if isfield(options, simOptions{i})
      simulator.(simOptions{i})  = options.(simOptions{i});
      options                    = rmfield(options, simOptions{i});
   end
end
%%
% =========================================================================
% Parallel computing management
% =========================================================================
% if isfield(options{1},'parallel')
%    %
%    if isa(options{1}.parallel,'logical') && options{1}.parallel
%       options{1}.parallel = struct('target', 'local');
%    end
%    %
%    if ~isfield(options{1}.parallel,'isRunning')
%       options{1}.parallel           = genetic.parallel.fillParallelOptions(options{1}.parallel);
%       options{1}.parallel.isRunning = true;
%       [xopt, fopt, info]            = genetic.parallel.executeParallel({fun, nobj}, xDim, methods, constraints, options);
%       bestCurrent.normalStop        = true;
%       return;
%    else
%       globalComm              = ompi.init('global');
%       comm                    = globalComm.create('evalComm', (1:globalComm.numtasks));
%       %
%       if comm.rank > 1
%          while true
%             params   = comm.receive_src(1);
%             cmd      = params{1};
%             if cmd == 0
%                ompi.finalize(0);
%             end
%             fun      = params{2};
%             xc       = params{3};
%             i        = params{4};
%             y        = fun(xc);
%             comm.send(1, {y, i});
%          end
%       end
%    end
%    %
%    options{1} = rmfield(options{1},'parallel');
% else
%    comm = [];
% end
%%
% =========================================================================
% Call to the optimization method(s)
% =========================================================================
%
%
optimizer            = genetic.optimizer.get(method, simulator, constraints, options);
%
optimizer            = optimizer.wrapper4constraints(constraints, options);
% 

[xopt, fopt, info]   = optimizer.apply(x0);
%
info.elapsedTime     = toc(tStart);
%
optimizer.printReport(info);
end