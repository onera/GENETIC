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

function paraConf = getLocalParallelOptions(userOptions)
paraConf            = struct('target','local','userName','','nCores',[],'isRunning',false);
% Local system identification
[userName, nCores]  = genetic.parallel.getLocalData();
paraConf.userName   = userName;
if nargin > 0 
    paraConf = genetic.tools.fillEmptyFields(userOptions,paraConf);
end
if isempty(paraConf.nCores)
    paraConf.nCores = nCores;
end
if paraConf.nCores > nCores 
    warning('Not enough cores locally, decreasing to %d',nCores)
    paraConf.nCores = nCores;
end

end 