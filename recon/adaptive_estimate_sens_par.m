% sens = adaptive_estimate_sens_par(varargin)
%
% Function to simplify running adaptive_estimate_sens() with parallel toolbox computing.
% Inputs are the same as adaptive_estimate_sens.m, except 'parts' and 'part' are taken care of.
%
% Written by Joseph G. Woods, QIS Lab, Bern, April 2025

function sens = adaptive_estimate_sens_par(varargin)

p = gcp;                 % Get the current parallel pool object
nWorkers = p.NumWorkers; % Find the number of workers available

% Extract inputs
p = inputParser;
p.addParameter('data',  []);
p.addParameter('kernel',5);
p.addParameter('thresh',0,  @isscalar);
p.addParameter('verbose',false, @islogical);
p.parse(varargin{:});
data    = p.Results.data;
kernel  = p.Results.kernel;
thresh  = p.Results.thresh;
verb    = p.Results.verbose;
clear p;

% Save size imformation and 
dims   = size(data);
N      = prod(dims(1:3));
ncoils = dims(4);
% data   = reshape(data, N, 1, 1, ncoils);
% 
% % Set up parallel job information and voxel indices for each worker to work on
% % and split up data to reduce memory overheads during parfor loop
% n   = ceil(N/nWorkers);
% idx      = cell(nWorkers,1);
% dataCell = cell(nWorkers,1);
% for ii = 1:nWorkers
%     idx{ii} = (ii-1)*n+1:min(ii*n, N);
%     dataCell{ii} = data(idx{ii},1,1,:);
% end
% clear data
% 
% % Run adaptive_estimate_sens, splitting it into nWorkers jobs
% sensCell = cell(nWorkers,1);
% maskCell = cell(nWorkers,1);
% parfor_progress(nWorkers); % Initialize
% parfor ii = 1:nWorkers
%     [sensCell{ii}, ~, maskCell{ii}] = ...
%         adaptive_estimate_sens('data', dataCell{ii}, 'kernel', kernel, 'thresh', thresh, 'verb', verb);
% end
% clear dataCell

% Run adaptive_estimate_sens, splitting it into nWorkers jobs
sensCell = cell(nWorkers,1);
idxCell  = cell(nWorkers,1);
maskCell = cell(nWorkers,1);
parfor_progress(N); % Initialize
for ii = 1:nWorkers
    [sensCell{ii}, idxCell{ii}, maskCell{ii}] = ...
        adaptive_estimate_sens('data', data, 'parts', nWorkers, 'part', ii, ...
                               'kernel', kernel, 'thresh', thresh, 'verbose', verb);
end
parfor_progress(0);

% Now sort it back into the shape used by adaptive_estimate_sens()
sens = zeros([ncoils, dims(1:3)]);
mask = zeros([1     , dims(1:3)]);
for ii = 1:nWorkers
    sens(:,idxCell{ii}) = sensCell{ii};
    mask(  idxCell{ii}) = maskCell{ii};
end
clear sensCell idxCell maskCell

% Threshold and reshape the sensitivity maps as done for nargin==1 in adaptive_estimate_sens.m
sens = permute(sens.*(abs(mask)>thresh*max(abs(mask(:)))),[2,3,4,1]);

end