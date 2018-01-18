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
function [pathToBench, EXT] = getPath()
geneticPath = genetic.tools.getPath();
basePath    = [geneticPath,filesep,'+bench',filesep,'+data',filesep];
pathToBench = {[basePath,'+mono',filesep,'+academic',filesep],...
               [basePath,'+mono',filesep,'+automatic',filesep],...
               [basePath,'+multi',filesep,'+academic',filesep],...
               [basePath,'+usr',filesep]};
EXT         = {'.m','.m','.m','.mat'};            
end