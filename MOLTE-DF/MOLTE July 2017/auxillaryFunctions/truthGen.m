function [f, beta_W,numD] = truthGen(problemClass, noiseRatio, varargin)
   if nargin == 2
       [f, beta_W,numD]=eval('caller',[problemClass, '(noiseRatio)']);
   else
        [f, beta_W,numD]=eval('caller',[problemClass, '(noiseRatio,varargin{1},varargin{2},varargin{3})']);
   end
end

