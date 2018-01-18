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
classdef mono <  matlab.unittest.TestCase
   properties
      meth = {'pso','cmaes','mdSearch','simplex','gss','linesearch'};
      n
      f
      x0
      f0
      cst
   end
   %
   methods(TestMethodSetup)
      function createObjectiveFunction(testCase)
         testCase.n  = 10;
         testCase.f  = @(x) norm(x)^2;
         testCase.x0 = 10*rand(testCase.n,1);
         testCase.f0 = testCase.f(testCase.x0);
         Aeq         = zeros(1,testCase.n);
         Aeq(1)      = 1;
         beq         = 1;
         testCase.cst = struct('Aeq',Aeq,'beq',beq);
      end
   end
   methods(Test)
      function allMethodsRun(testCase)
         errors         = [];
         opt.maxFunEval = 200;
         for i = 1:length(testCase.meth)
            m = testCase.meth{i};
            try
               [xopt, fopt] = genetic.min(testCase.f, testCase.n, m, opt);
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
      %%
      function allWorkingMethodsDoNotDegrade(testCase)
         opt.maxFunEval = 1000;
         errors = [];
         for i = 1:length(testCase.meth)
            m = testCase.meth{i};
            try
               [xopt, fopt]   = genetic.min(testCase.f, testCase.x0, m, opt);
               if fopt > testCase.f0
                  if isempty(errors)
                     sep = '';
                  else
                     sep = ',';
                  end
                  errors = [errors,sep,m];
               end
            catch
            end
         end
         testCase.verifyEmpty(errors,['Method not improving initial guess:',errors]);
      end
   end
end