function paretoFront(varargin)
%
Y        = {}; 
nobj     = [];
nFronts  = 0;
toColor = [];
for i = 1:length(varargin)
   Yi = varargin{i};
   if isa(Yi, 'double')
      nFronts        = nFronts + 1;
      if isempty(nobj)
         nobj        = size(Yi,1);
         if nobj > 3
            error('Displaying Pareto front is only supported for dimension 2 and 3')
         end
      else
         if size(Yi,1) ~= nobj
            error('Incoherent number of objectives with front %d',nFronts)
         end
      end
      for j = 1:nobj
         Y{end+1} = Yi(j,:);
      end
   end
   
   if i < length(varargin)
      if isa(varargin{i+1},'char')
         Y{end + 1}  = varargin{i+1};
         i           = i+1;
      else
         toColor = [toColor,i];
      end
   end
end
% marks = {'.','o','x','+','*','s','d','v','^','<','>','p','h'};
% ctr = 0;
% for i = sort(toColor,'descend')
%    
% end
%
if nobj == 2
   plot(Y{:})
else
   plot3(Y{:})
end

   
end