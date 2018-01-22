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
function [constraints, options] = separateConstraintsAndOptions(args)
constraints = [];
options     = [];
if length(args) == 1
   % if there is only one additional input argument, it may be the
   % constraints structure or the options structure.
   % To distinguish the two cases, the fields of the argument are tested.
   constraintsFields = {'A','b','Aeq','beq','c','ceq','xMin','xMax'};
   fields            = fieldnames(args{1});
   if ~isempty(intersect(fields,constraintsFields)) && isempty(setdiff(fields, constraintsFields))
      % If the argument contains at least one of the field specific to
      % constraints and no additional field, then it is the constraints
      % structure
      constraints = args{1};
   else
      % otherwise, it is the options structure
      options     = args{1};
   end
elseif length(args) == 2
   % If there are two arguments, then the options structure is the last
   % one.
   constraints = args{1};
   options     = args{2};
elseif length(args)> 2
   error('Too many input arguments')
end
end