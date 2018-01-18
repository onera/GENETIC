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
function out = fillEmptyFields(in, default)
% This routine adds and fills the fields that are missing in the structure 
% "in" in comparison to the structure "default".
% example :
% s = struct('a',1,'b',2);
% ds = struct('a',-1,'b',-2,'c',-3);
% genetic_fillEmptyFields(s,ds)
% returns struct('a',1,'b',2,'c',-3);
if isempty(in)
   in = struct();
end
out            = in;
fieldNames     = fieldnames(default);
fieldsToFill   = find(~isfield(in, fieldNames));
for i = 1:length(fieldsToFill)
   idField = fieldsToFill(i);
   out.(fieldNames{idField}) = default.(fieldNames{idField});
end

end