function hv = hyperVolume(Y, ref, method, varargin)
% Example
% % In 2 dim
% n      = 200;
% theta  = linspace(-3*pi/2,-pi,n)';
% PF     = [1+cos(theta) 1-sin(theta)]';
% hv     = genetic.population.metrics.hyperVolume(PF);
% % hv should tend to pi/4

idx = genetic.population.mgroup.getNonDominatedIndividuals(Y);
Y   = Y(:,idx);

if nargin < 2 || isempty(ref)
   ref = max(Y,[],2);
end
%
if nargin < 3
    if size(Y,1) == 2
        method = 'hso';
    else
        method = 'leb';
    end
end
%
Y           = unique(Y','rows')';
[nObj, n]   = size(Y);
%
hv = [];
switch method      
   case 'hso'
      if nObj == 2
         hv = hv2D(Y, ref);
         return
      end
      if nObj == 3
%          hv = hv3D(Y, ref);
         return
      end
   case {'leb','lebesgue'}
      hv = hvLeb(Y, ref);
   case {'mc','monte-carlo','approx'}
      hv = hvEst(Y, ref, varargin{:});
   otherwise
      error('not implemented')
end
   

end


function hv = hv2D(Y, ref)
n  = size(Y,2);
% Sort Y based on the first objective in ascending order
Y  = sortAccordingToDim(Y, 1);
hv = sum((ref(1) - Y(1,:)) .* ([ref(2), Y(2,1:n-1)] - Y(2,:)));
end

% function hv = hv3D(Y, ref)
% Paquete's algorithm
% [nObj, n]   = size(Y);
% Sort according to third dimension
% Y           = sortAccordingToDim(Y, 3);
% %
% hv          = 0;
% p           = Y(:,1);
% ps          = Y(:,end);
% %
% area        = p(1) * p(2);
% z           = p;
% for i =  
% end
% end


function Ys = sortAccordingToDim(Y, d)
Ys = sortrows(Y',d)';
end

function hv = hvLeb(L, ref)
% M. Fleischer - The measure of Pareto optima: Applications to
% multi-objective metaheuristics. Technical report, 2002.
[nObj, n]   = size(L);
% L           = sortAccordingToDim(Y,1);
hv          = 0;
while n > 1
   p1          = L(:,1);
   sv          = repmat(p1, 1, nObj); % spawned vectors
   delta       = bsxfun(@gt, L , p1);
   b           = zeros(nObj,1);
   for i = 1:nObj
      if any(delta(i,:))
         b(i)  = min(L(i,delta(i,:)));
      else
         b(i) = ref(i);
      end
      sv(i,i)  = b(i);
   end
   lopOffVol   = prod(abs(b - p1));
   hv          = hv + lopOffVol;
   L           = ndFilter(L(:,2:end), sv, ref);
   % Update number of points
   n           = size(L,2);
end
% % Last iteration
% hv = hv + prod(ref - L(:,1));


end

function L = ndFilter(L, sv, ref)
addToL   = false(1,size(sv,2));
for i = 1:size(sv,2)
   % for each spawned vector
   svi = sv(:,i);
   if all(svi ~= ref)
      addToL(i) = true;       
       svIsDominated  = all(bsxfun(@le, L, svi),1) & any(bsxfun(@ne,L,svi),1);
       if sum(svIsDominated) == 0
          addToL(i) = true;
      end
   end
end
L = [sv(:,addToL),L];
end



function hv = hvEst(Y, ref, N)
   % Estimation of the hypervolume based on Monte-Carlo
   if nargin < 3
      N = 1000;
   end
   [nObj,n] = size(Y);
   lb       = min(Y, [], 2);
   % Normalization of Y
   Y        = bsxfun(@times,Y,1./ref);
   %
   C           = rand(nObj, N);
   fDominated  = false(N,1);
   fcheck      = all(bsxfun(@gt, C, lb), 1);
   %
   for k = 1:n
      if any(fcheck)
         f                    = all(bsxfun(@gt, C(:,fcheck), Y(:,k)), 1);
         fDominated(fcheck)   = f;
         fcheck(fcheck)       = ~f;
      end
   end
   %
   hv          = sum(fDominated)/N;
end