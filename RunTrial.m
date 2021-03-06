function [data, timedOut, quit] = RunTrial(data, trial, window, varargin)

% To be run at the beginning of every trial
%
% data   -  for stim parameters and saving responses
% trial  -  trial num (0 -> practice)
% window -  for showing stimuli

%% Prepare Trial

% set cue text size
Screen('TextSize', window, data.exp.cueTextSize);

% Set random seed
randSeed = rng('shuffle');

% Select Parameters
c = binornd(1, 0.5) + 1; % contrast
contrast = data.stimuli.contrastVal(c);

v = binornd(1, 0.5) + 1; % variability
variability = data.stimuli.variabilityVal(v);

if binornd(1, 0.5) % Generative Distribution
    genDist = 3;
else
    genDist = -3;
end

genMean = randn*data.stimuli.genStd + genDist % orientation mean
correct = round((sign(genMean)+1)/2); % correct -> sign of orientation mean

orientations = NaN(data.stimuli.nSamples, 1);
for i=1:data.stimuli.nSamples % orientations for each gabor
        orientations(i) = randn*variability - genMean;
end

%orientations = orientations + (genMean - mean(orientations));

% Generate Cue
if trial
    if data.response.isCuedBlock(:,trial)==0 % uncued block
        cue = 0;
        cueLet = 'N';
    else                                     % cued block
        cue = sign(binornd(1, abs(correct-0.25)) - 0.5);
        switch cue
            case -1
                cueLet = 'L';
            case 1
                cueLet = 'R';
        end
    end
else % practice trials
    cue = binornd(2, 0.5) - 1; % make this probabilistic as well
    switch cue
        case -1
            cueLet = 'L';
        case 0
            cueLet = 'N';
        case 1
            cueLet = 'R';
    end
end

% Generate Stimuli
gabors = ones(data.stimuli.gaborDimPix, data.stimuli.gaborDimPix, data.stimuli.nSamples);

for i=1:data.stimuli.nSamples % Generate Gabors
   gabors(:,:,i) = 0.5+SCreateGaborINT(data.stimuli.gaborDimPix, data.stimuli.envelopDev,...
       orientations(i), data.stimuli.frequency, data.stimuli.phase(i, max(trial, 1)),...
       contrast, data.stimuli.noiseAmp) .* data.stimuli.aperture;
end

%% Show trial

% show fixation point
Screen('FillOval', window, [0 0 0], data.stimuli.fixBox);
vbl = Screen('Flip', window); % get time

% show cue
DrawFormattedText(window, cueLet, 'center', 'center', [0 0 0]);
vbl = Screen('Flip', window, vbl + data.exp.fixShowTime, 1); % don't clear screen

% show Gabors
for i=1:data.stimuli.nSamples
   Screen('PutImage', window, gabors(:,:,i), data.stimuli.gratingBoxes(i,:)');
end

vbl = Screen('Flip', window, vbl + data.exp.cueShowTime); % show stimuli

DrawFormattedText(window, cueLet, 'center', 'center', [0 0 0]); % keep cue during response
vbl = Screen('Flip', window, vbl + data.exp.stimShowTime); % end stimuli

%% Reports

% Get Response
allowResponse = 1; % loop until response
startSec = GetSecs; % get time
quit = 0; % default no quit

while allowResponse
    [kbDown, sec, kbKey] = KbCheck; % get kb info and time

    % handle timeout
    if (sec-startSec) >= data.exp.responseTimout
        allowResponse = 0;
        timedOut = 1;
    end
    
    if kbDown==1
        if kbKey(data.exp.leftKey)==1 % CCW
            allowResponse = 0;
            timedOut = 0;
            response = 0
            RT = sec - startSec;
        elseif kbKey(data.exp.rightKey)==1 % CW
            allowResponse = 0;
            timedOut = 0;
            response = 1
            RT = sec - startSec;
        elseif kbKey(data.exp.escapeKey)==1 % Quit
            allowResponse = 0;
            timedOut = 0;
            quit = 1;
            
            return % exit func
 
        end % kbKey        
    end % kbDown
end % while allowResponse


% set normal font size
Screen('TextSize', window, data.exp.directionTextSize);


% Get Confidence Report
allowResponse = 1; % loop until response
quit = 0; % default no quit

while allowResponse && ~timedOut
    [kbDown, sec, kbKey] = KbCheck; % get kb info and time
    
    DrawFormattedText(window, 'How Confident Are You?\n8 - High Confidence\n5 - Medium Confidence\n2 - Low confidence', 'center', 'center', [0 0 0]);
    vbl = Screen('Flip', window); 
    
    if kbDown==1  
        if kbKey(data.exp.high)==1 % High Confidence
            allowResponse = 0;
            confidence = 1;
        elseif kbKey(data.exp.medium)==1 % Med Confidence
            allowResponse = 0;
            confidence = 0;
        elseif kbKey(data.exp.low)==1 % Low Confidence
            allowResponse = 0;
            confidence = -1;
        elseif kbKey(data.exp.escapeKey)==1 % Quit
            allowResponse = 0;
            quit = 1;
            return
        end % kbKey        
    end % kbDown
end% while allowResponse


%% Feedback and Sava Data

if timedOut
    DrawFormattedText(window, sprintf('Please answer within %i seconds', data.exp.responseTimout),...
        'center', 'center', [1 1 1]);
else % successful trial
    accuracy = correct==response;
    
    switch accuracy % generate feedback text
        case 0
            DrawFormattedText(window, 'Incorrect!', 'center', 'center', [1 0 0]);
        case 1
            DrawFormattedText(window, 'Correct!', 'center', 'center', [0 1 0]);
    end
end

vbl = Screen('Flip', window); % show feedback
Screen('Flip', window, vbl + data.exp.feedbackShowTime);


% Save Data
if trial && ~timedOut
    data.response.randSeed(:, trial) = randSeed.Seed;
    
    data.response.contrast(:,trial) = contrast;
    data.response.variance(:,trial) = variability;
    data.response.correct(:,trial) = correct;
    
    data.response.orientationMean(:,trial) = genMean;
    data.response.trueOrientaions(:,trial) = orientations;
    data.response.cue(:,trial) = cue;
    
    data.response.responseRight(:,trial) = response;
    data.response.accuracy(:,trial) = accuracy;
    data.response.reactionTime(:,trial) = RT;  
  	data.response.confidence(:,trial) = confidence;    
end % end save data

end % end RunTrial func


%% Functions from Target Paper

function [patch] = SCreateGaborINT(siz,envelopedev,angle,frequency,phase,contrast,noiseamp)
    if nargin < 6, contrast = 1; end
    if nargin < 5, error('Not enough input arguments.'); end
    
    [xint,yint] = meshgrid((1:siz)-(siz/2));%-(siz/2));%xint and yint are the internal versions
    
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