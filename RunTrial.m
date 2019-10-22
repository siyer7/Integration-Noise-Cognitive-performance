function [data] = RunTrial(data, trial, varargin)

% To be run at the beginning of every trial
%
% Outline
% 1. Generate stimulus and cue (save)
% 2. Show stimulus
% 3. Response
% 4. Feedback
% 5. Save and return data

% Select Parameters
c = binornd(1, 0.5) + 1;
contrast = data.stimuli.contrastVal(c);

v = binornd(2, 0.5) + 1;
variability = data.stimuli.variabilityVal(v); % actuall std dev

if binornd(1, 0.5)
    correct = 1;
    genDist = 3;
else
    correct = 0;
    genDist = -3;
end

genMean = randn*data.stimuli.genStd + genDist;

orientations = NaN(data.stimuli.nSamples, 1);
for i=1:data.stimuli.nSamples
    orientations(i) = randn*variability + genMean;
end

if binornd(1, 0.5)
    cue = 0; % uninformative cue
else
    cue = sign(binornd(1, abs(correct-0.25)) - 0.5);
end

% Generate and Show Stimuli
gabors = ones(data.stimuli.gaborDimPix, data.stimuli.gaborDimPix, data.stimuli.nSamples);
for i=1:data.stimuli.nSamples
   gabors(:,:,i) = SCreateGaborINT(data.stimuli.gaborDimPix, data.stimuli.envelopDev,...
       orientations(i), data.stimuli.frequency, data.stimuli.phase(i, max(trial, 1)),...
       contrast, data.stimuli.noiseAmp) .* data.stimuli.aperture;
   
   gabors(:,:,i) = gabors(:,:,i) + min(gabors(:,:,i));
   indx = Screen('MakeTexture', window, gabors(:,:,i));
   Screen('PutImage', window, indx, data.stimuli.gratingBoxes(i,:)');
end

Screen('Flip', window);

KbStrokeWait;

% Get Response

% Save Data
if trial
    data.response.contrast(:,trial) = contrast;
    data.response.variance(:,trial) = variability;
    data.response.correct(:,trial) = correct;
    
    data.response.orientationMean(:,trial) = genMean;
    data.response.trueOrientaions(:,trial) = orientations;
    data.response.cue(:,trial) = cue;
    
    %data.response.responseRight(:,trial) = 0;
    %data.response.accuracy(:,trial) = NaN(1, data.exp.numTrials);
    %data.response.reactionTime(:,trial) = NaN(1, data.exp.numTrials);  
    
end


end

%-----------------------------
% Functions from Target Paper
%-----------------------------

function [patch] = SCreateGaborINT(siz,envelopedev,angle,frequency,phase,contrast,noiseamp)
    if nargin < 6, contrast = 1; end
    if nargin < 5, error('Not enough input arguments.'); end
    
    [xint,yint] = meshgrid((1:siz));%-(siz/2));%xint and yint are the internal versions
    
    patch = 0.5*contrast*cos(2*pi*(frequency*(sin(pi/180*angle)*xint+cos(pi/180*angle)*yint)+phase));

    patch = patch + CreateSmoothedNoiseINT(siz,frequency,noiseamp);
    
    envel = NormalPDF(sqrt(xint.^2+yint.^2),0,envelopedev);
    
    patch = patch.*(envel./(max(envel(:))));
end

function [patch] = CreateSmoothedNoiseINT(siz,frequency,dev)
    kerneldev = 1/frequency;
    
    if nargin < 3, dev = 1; end
    if nargin < 2, error('Not enough input arguments.'); end
    
    kernelsiz = ceil(kerneldev);
    kernellim = 0.5*kernelsiz/kerneldev;
    
    kernel = NormalPDF(linspace(-kernellim,+kernellim,kernelsiz),0,1);
    kernel = kernel'*kernel;
    
    patch = conv2(randn(siz),kernel,'same');
    patch = (patch/sqrt(sum(kernel(:).^2)))*dev;
end