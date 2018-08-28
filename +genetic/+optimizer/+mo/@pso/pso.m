% PSO - Particle Swarm Optimisation
%
%  PSO is an evolutionary optimization scheme originaly inspired by the
%  movement of bird flock or fish school.
%
%  A population (swarm) of candidate solutions (particules) moves in the
%  search-space according to some simple formulae. These formulae
%  involve a randomized combination of various elementary movements: 
%  a movement induced by the inertia of the particle, a movement towards
%  the best particle in the swarm, a movement towards the best 
%  neighbour, etc.
%
%  The neighbours are determined according to some network topology in 
%  the swarm (Von Neuman, full, ring).
%
%  For further informations, see e.g. [1] and references therein.
%
% Parameters
%  topologyType: type of topology to connect the particules
%                'vonn' (default), 'full', 'ring'
%
%  randDistrib : type of distribution for the velocity update
%                'rectangular' (default), 'spheric', 'gaussian',
%                'local1', 'local2', 'pivot1', 'pivot2', 'pivot3', 
%                'pivot4', 'pivot5', 'pivot6'
%
% References
%  [1] M. Clerc. Particle swarm optimization. ISTE. 2006.
%

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
classdef pso < genetic.optimizer.mono & genetic.optimizer.simpleScheme
   properties
      % Default options
      topologyType      = 'vonn';
      randDistrib       = 'rectangular';
      phiConstant       = 2.13;
      %
   end
   properties (Access = private)
      consFactor
      cMax
      speedField
      calls
   end
   methods
      %% constructor
      function self = pso(varargin)
         self@genetic.optimizer.mono(varargin{:});
         self@genetic.optimizer.simpleScheme();
         %
         self.methodName         = 'pso';
         self.longName           = 'Particle Swarm optimisation';
         self.needFiniteBounds   = true;
         %
         self.calls              = {@(I) I.value, @(I) I.best.value, @(I) I.getBestNeighbour().value, @(I) I.best.objective, @(I) I.getBestNeighbour().best.objective};
         self.consFactor         = 1/(self.phiConstant-1+sqrt(self.phiConstant^2-2*self.phiConstant));
         self.cMax               = self.consFactor * self.phiConstant;
         self.speedField         = repmat(self.xMin,1, self.popSize) + rand(self.xDim, self.popSize).*repmat((self.xMax - self.xMin)/2, 1, self.popSize );
      end
      %% evolve
      % Return a new group by applying specific mechanisms (proper to this optimizer)
      function group = evolve(obj, group)
         xDim     = group.xDim;
         popSize  = group.len;
         % velocity vector(s) field computation
         data     = group.applyToIndividuals(obj.calls);
         % access to current particles's positions
         x        = data{1};
         pBest    = data{2};
         nBest    = data{3};
         valPBest = repmat(data{4},xDim,1);
         valNBest = repmat(data{5},xDim,1);
         % random distribution construction
         switch obj.randDistrib
            case 'rectangular'
               randPBest      = obj.cMax*rand(xDim,popSize);
               randNBest      = obj.cMax*rand(xDim,popSize);
               % update velocity vectors field
               obj.speedField = obj.consFactor*obj.speedField + randPBest.*(pBest - x) + randNBest.*(nBest - x);
               % update particles' positions
               newX           = x + obj.speedField;
            case 'spheric'
               % formula sphere volume in dimension d: Vd = (pi^(d/2)*radius^d)/gamma(1+d/2)
               radius    = sqrt(3)*obj.cMax/2;
               tmp1      = randn(popSize,xDim);
               tmp2      = sum(tmp1.^2,2);
               X         = tmp1.*repmat(radius*(gammainc(tmp2/2,xDim/2).^(1/xDim))./sqrt(tmp2),1,xDim);
               randPBest = obj.cMax/2 + X';
               tmp1      = randn(popSize,xDim);
               tmp2      = sum(tmp1.^2,2);
               X         = tmp1.*repmat(radius*(gammainc(tmp2/2,xDim/2).^(1/xDim))./sqrt(tmp2),1,xDim);
               randNBest = obj.cMax/2 + X';
               % update velocity vectors field
               obj.speedField = obj.consFactor*obj.speedField + randPBest.*(pBest - x) + randNBest.*(nBest - x);
               % update particles' positions
               newX     = x + obj.speedField;
            case 'gaussian'
               randPBest     = abs(obj.cMax/2 + obj.cMax*randn(xDim,popSize)/6);
               randNBest     = abs(obj.cMax/2 + obj.cMax*randn(xDim,popSize)/6);
               % update velocity vectors field
               obj.speedField = obj.consFactor*obj.speedField + randPBest.*(pBest - x) + randNBest.*(nBest - x);
               % update particles' positions
               newX     = x + obj.speedField;
            case 'local1'
               % update particles' positions
               newX     = nBest + abs(nBest - x).*randn(xDim,popSize)/2;
            case 'local2'
               % update particles' positions
               newX     = (nBest - x)/2 + abs(nBest - x).*randn(xDim,popSize)/2;
            case {'pivot1', 'pivot2', 'pivot3','pivot4','pivot5','pivot6'}
               pivotNumber = str2double(obj.randDistrib(6));
               % radii evaluation between <pBest> and <nBest>
               radii     = repmat(sqrt(sum((pBest - nBest).^2,1)),xDim,1);
               if all(pivotNumber ~= [4,5,6])
                  %~strcmp(obj.randDistrib,'pivot4')&&~strcmp(obj.randDistrib,'pivot5')&&~strcmp(obj.randDistrib,'pivot6')
                  tmp1      = randn(popSize,xDim);
                  tmp2      = sum(tmp1.^2,2);
                  X         = tmp1.*radii'.*repmat((gammainc(tmp2/2,xDim/2).^(1/xDim))./sqrt(tmp2),1,xDim);
               else
                  X         = radii'.*randn(popSize,xDim)/2.326/3;
               end
               randPBest = pBest + X';
               if all(pivotNumber ~= [4,5,6])
                  %~strcmp(obj.randDistrib,'pivot4')&&~strcmp(obj.randDistrib,'pivot5')&&~strcmp(obj.randDistrib,'pivot6')
                  tmp1      = randn(popSize,xDim);
                  tmp2      = sum(tmp1.^2,2);
                  X         = tmp1.*radii'.*repmat((gammainc(tmp2/2,xDim/2).^(1/xDim))./sqrt(tmp2),1,xDim);
               else
                  X         = radii'.*randn(popSize,xDim)/2.326/3;
               end
               randNBest = nBest + X';
               % update particle's position
               newX = (valPBest.*randPBest + valNBest.*randNBest)./(valPBest + valNBest);
               % additive gaussian noise
               if all(pivotNumber ~= [1,4])
                  %~strcmp(obj.randDistrib,'pivot1')&&~strcmp(obj.randDistrib,'pivot4')
                  if any(pivotNumber == [2,5])
                     %strcmp(obj.randDistrib,'pivot2')||strcmp(obj.randDistrib,'pivot5')
                     % per component
                     tmp = randn(xDim,popSize);
                  elseif any(pivotNumber == [3,6])
                     %strcmp(obj.randDistrib,'pivot3')||strcmp(obj.randDistrib,'pivot6')
                     % per individual
                     tmp = repmat(randn(1,popSize),xDim,1);
                  end
                  % noise matrix computation
                  noise = (valPBest - valNBest).*tmp./(valPBest + valNBest);
                  % newX construction
                  newX = (1 + noise).*newX;
               end
            %case 'ellipsoid'
               % TODO ??
%                else
%                    error('unknown
         end
         group.moveTo(newX);
      end
      
   end
end

