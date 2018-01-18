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
function x = voronoi(popSize, xMin, xMax, xDim, varargin1, varargin2)
if nargin == 4
   varargin1 =   [];
   varargin2 = 2000;
end
% Voronoi parameters definition <nIter> + <n> + <lambda> + <refError>
nIter    = 50;
n        = varargin2;
lambda   = 0.75;
refError = 0.45e-02*norm(abs(xMax-xMin));
% MATLAB 'if' loop on <varargin1>
if isempty(varargin1)
   x = genetic.population.init.dispatchIndividuals(popSize, xMin, xMax, rand(xDim,popSize));
else
   x = varargin1;
end
% MATLAB 'for' loop on integer <nIter>
for i = 1:nIter
   % Random generation of <n> samples with dimension <xDim>
   ptc = repmat(xMin,1,n) + repmat(abs(xMax-xMin),1,n).*rand(xDim,n);
   % Parameters <dist> + <nbwin> + <subset> initializations
   dist     = zeros(popSize,1);
   nbwin    = zeros(popSize,1);
   subset   = zeros(popSize,0);
   % MATLAB 'for' loop on integer parameter <n>
   for j = 1:n
      % MATLAB 'for' loop on <popSize>
      for k = 1:popSize
         dist(k) = norm(x(:,k)-ptc(:,j));
      end
      [~,l]                = min(dist);
      nbwin(l)             = nbwin(l) + 1;
      subset(l,nbwin(l))   = j;
   end
   % Variables <outFlag> + <epsMax> initializations
   outFlag  = 1;
   epsMax   = 0;
   % MATLAB 'for' loop <popSize>
   for j = 1:popSize
      % MATLAB 'if' loop on <nbwin(j)>
      if nbwin(j) ~= 0
         center = mean(ptc(:,subset(j,1:nbwin(j))),2);
         eps = center - x(:,j);
         r = norm(eps);
         % MATLAB 'if' loop on <r>
         if r > refError
            outFlag = 0;
         end
         % MATLAB 'if' loop on <r>
         if r > epsMax
            epsMax = r;
         end
         x(:,j) = x(:,j) + lambda*eps;
      end
   end
   % Breaking condition on <outFlag>
   if outFlag
      break
   end
end

end