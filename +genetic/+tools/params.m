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
function out = params(varargin)
persistent params
if isempty(params)
   params = struct('verbosity',0,'dbg',false);
end
if length(varargin) > 1 && mod(length(varargin),2) == 0
   in       = genetic.tools.keyValuePairs(varargin{:});
   params   = genetic.tools.fillStructure(params, in);
elseif length(varargin) == 1 && isa(varargin{1},'char')
   out      = params.(varargin{1});
end
end