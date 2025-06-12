% sens = adaptive_estimate_sens_par(varargin)
%
% Function to simplify running adaptive_estimate_sens() with parallel toolbox computing.
% Inputs are the same as adaptive_estimate_sens.m, except 'parts' and 'part' are taken care of.
%
% Written by Joseph G. Woods, QIS Lab, Bern, April 2025

function sens = adaptive_estimate_sens_par(varargin)

% Extract inputs
p = inputParser;
p.addParameter('data',  []);
p.addParameter('kernel',5);
p.addParameter('thresh',0,  @isscalar);
p.addParameter('verbose',false, @islogical);
%p.addParameter('maxNWorkers',inf, @isscalar);
p.parse(varargin{:});
data    = p.Results.data;
kernel  = p.Results.kernel;
thresh  = p.Results.thresh;
verb    = p.Results.verbose;
%maxNWorkers = p.Results.maxNWorkers;
clear p;

% pool = gcp;                 % Get the current parallel pool object
% nWorkers = pool.NumWorkers; % Find the number of workers available
% if nWorkers > maxNWorkers
%     nWorkers = maxNWorkers;
% end

% Save size imformation and 
dims   = size(data);
%N      = prod(dims(1:3));
ncoils = dims(4);

% Run adaptive_estimate_sens, splitting it into multiple jobs
% I use 4*nWorkers so we get an update roughly every 25% chunk completed
N = 100;
sensCell = cell(N,1);
idxCell  = cell(N,1);
maskCell = cell(N,1);
parfor_progress(N); % Initialize
parfor ii = 1:N
    [sensCell{ii}, idxCell{ii}, maskCell{ii}] = ...
        adaptive_estimate_sens('data', data, 'parts', N, 'part', ii, ...
                               'kernel', kernel, 'thresh', thresh, 'verbose', verb);
    parfor_progress;
end
parfor_progress(0);

% Now sort it back into the shape used by adaptive_estimate_sens()
sens = zeros([ncoils, dims(1:3)]);
mask = zeros([1     , dims(1:3)]);
for ii = 1:N
    sens(:,idxCell{ii}) = sensCell{ii};
    mask(  idxCell{ii}) = maskCell{ii};
end
clear sensCell idxCell maskCell

% Threshold and reshape the sensitivity maps as done for nargin==1 in adaptive_estimate_sens.m
sens = permute(sens.*(abs(mask)>thresh*max(abs(mask(:)))),[2,3,4,1]);

end