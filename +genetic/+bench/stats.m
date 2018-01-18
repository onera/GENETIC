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
classdef stats < handle
   properties
      %
      n
      %
      fopt
      elapsedTime
      nEval
      nImprove
      maxCstViolation
      %
      nSuccess    = 0;
      nError      = 0;
      successRate = 0;
      errorRate   = 0;
      %
   end
   methods
      function self = stats(res)
         self.fopt            = genetic.bench.stats.empty();
         self.elapsedTime     = genetic.bench.stats.empty();
         self.nEval           = genetic.bench.stats.empty();
         self.nImprove        = genetic.bench.stats.empty();
         self.maxCstViolation = genetic.bench.stats.empty();
         self.n               = length(res);
         for i = 1:self.n
            self.add(res{i},i);
         end
         self.successRate     = self.nSuccess / self.n;
         self.errorRate       = self.nError / self.n;
      end
      %
      function add(self, r, i)
         %
         minMaxFields = {'fopt','elapsedTime','nEval','nImprove','maxCstViolation'};
         for j = 1:length(minMaxFields)
            f                 = minMaxFields{j};
            self.(f).min      = min([self.(f).min, r.(f)]);
            self.(f).max      = max([self.(f).max, r.(f)]);
            if i >  1
               self.(f).var   = (i-2)/(i-1) * self.(f).var + (r.(f) - self.(f).mean)^2/i;
            end
            self.(f).std      = sqrt(self.(f).var);
            self.(f).mean     = self.(f).mean + 1/i * (r.(f) - self.(f).mean); % Iterative computation of the mean value
         end
         %
         if r.success
            self.nSuccess     = self.nSuccess + 1;
         end
         %
         if r.error
            self.nError       = self.nError + 1;
         end
      end
   end
   %
   methods(Static)
      function out = empty()
         out = struct('min',inf,'max',-inf,'mean',0,'std',[],'var',0);
      end
   end
end
   