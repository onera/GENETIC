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
function list()
% GENETIC.BENCH.LIST lists the existing benchmarks.
[pathToBench,EXT] = genetic.bench.getPath();
head              = {'Mono objective (academic)','Mono objective (automatic)','Multi objective (academic)', 'User defined'};
opt.alignment     = {'l','l','l','l','l','l'};
ncols             = 6;
for i = 1:length(pathToBench)
   tmp      = dir([pathToBench{i},'*',EXT{i}]);
   doPrint = ~isempty(tmp);
   if doPrint
      nrows    = max(ceil(length(tmp)/ncols),1);
      entries  = cell(nrows,ncols);
      n        = 0;
      for j = 1:ncols
         for k = 1:nrows
            n = n + 1;
            try
               post           = '';
               benchName      = strrep(tmp(n).name,EXT{i},'');
               b              = genetic.bench.load(benchName, [], true);                  
               if ~isempty(b)
                  post        = ['x:',num2str(b.xDim)];
               else
                  b           = genetic.bench.load(benchName, 4, true);
               end
               if b.fDim > 1
                  if ~isempty(post)
                     post = [post,','];
                  end
                  post  = [post, 'y:',num2str(b.fDim)];
               end
               if ~isempty(post)
                  post = [' (',post,')'];
               end
               entries{k}{j}  = [benchName,post];
            catch
               entries{k}{j} = ' ';
            end
         end
      end
      if nrows == 1
         entries = {entries{1}}; % the cell is weirdly furnished when nrows == 1
      end
      tab = genetic.tools.tabular(head{i},[], entries, opt);
      fprintf(tab)
   end
end
end