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
classdef multi <  matlab.unittest.TestCase
   properties
      meth = {'mombi2','mopso','nsga2','nsga3','pesa2','rveastar','spea2','tdea'};
      n
      f
      f0
      nobj
   end
   methods(TestMethodSetup)
      function createObjectiveFunction(testCase)
         testCase.n     = 10;
         A              = rand(testCase.n,testCase.n);
         b              = rand(testCase.n,1);
         testCase.f     = @(x) [norm(A*x-b)^2;norm(x,inf)];
         testCase.nobj  = 2;
      end
   end
   methods (Test)
      function allMethodsRun(testCase)
         errors         = [];
         opt.maxFunEval = 200;
         for i = 1:length(testCase.meth)
            m = testCase.meth{i};
            try
               [xopt, fopt] = genetic.min({testCase.f,testCase.nobj}, testCase.n, m, opt);
            catch ME
               if isempty(errors)
                  sep = '';
               else
                  sep = ',';
               end
               errors = [errors,sep,m];
            end
         end
         testCase.verifyEmpty(errors,['Problem with methods:',errors]);
      end      
   end
end