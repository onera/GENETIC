% GENETIC.MAX - Main interface for maximisation.
% It can be used as <a href="matlab: help genetic.min">genetic.min</a>.

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
function [xopt, fopt, info] = max(fun, xDim, method, varargin)
% Extracting some options to adapt them for minimisation
[constraints, options]  = genetic.tools.separateConstraintsAndOptions(varargin);
% New objective function is the opposite of the initial one
newFun = @(x) -fun(x);
% If the gradient if provided, then the wrapper must concern also the
% gradient output
if isfield(options,'gradObj') && options.gradObj
   newFun = @(x) wrapf_(x, fun);
end
% If a target in objective is provided, its sign must be modified
if isfield(options, 'targetY')
   options.targetY = -options.targetY;
end

genetic.tools.params('ffactor', -1);
%
[xopt, fopt, info]      = genetic.min(newFun, xDim, method, constraints, options);
fopt                    = -fopt;
end

function [y, g] = wrapf_(x, f)
if nargout == 1
   y     = -f(x);
else
   [y,g] = f(x);
   y     = -y;
   g     = -g;
end
end
