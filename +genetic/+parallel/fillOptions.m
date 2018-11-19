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

function opt = fillOptions(opt)

if isa(opt,'logical')
    % parallel option is just a boolean...
    if opt
        % if true then the target is local
        opt     = struct('');
        target  = 'local';
    else
        % if false, then it is not parallelised and one returns
        opt = [];
        return
    end
elseif isa(opt, 'struct')
    % parallel option is a structure of options
    target = opt.target;
end
    

switch target
   case 'local'
    opt = genetic.parallel.getLocalParallelOptions(opt);
end
%

end