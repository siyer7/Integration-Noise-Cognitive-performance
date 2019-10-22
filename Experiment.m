function [success, data] = Experiment(varargin)

% Main Code
%
% Outline
% 
% 1. Get subj data and parameters
% 2. Prep Psychtoolbox
% 3. Run Experiment
% 4. Save Data

% Clear Workspace
close all;
sca;
clear vars;

% Get Subject Data and Parameters
data = GetData();

directory = fullfile(pwd);
fileName = fullfile(directory, [data.participant.name '.mat']);


% Prepare Psychtoolbox Stuff
PsychDefaultSetup(2);

screenNum = max(Screen('Screens'));

black = [ 0  0  0];
grey =  [.5 .5 .5];
white = [ 1  1  1];

[window, windowRect] = PsychImaging('OpenWindow', screenNum, grey, [310, 140, 1610, 940], 32, 2,...
        [], [],  kPsychNeed32BPCFloat);

[xCenter, yCenter] = RectCenter(windowRect);
[xWidth, yHeight] = Screen('WindowSize', window); 
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

ifi = Screen('GetFlipInterval', window);

topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel); % ERROR idk why :(

Screen('TextSize', window, data.exp.directionTextSize);

% Prepare Stimulus Stuff
data.stimuli.aperture = CreateCircularApertureINT(data.stimuli.gaborDimPix, 1);

posX = data.stimuli.posX';
posY = data.stimuli.posY';
gaborDim = data.stimuli.gaborDimPix;
data.stimuli.gratingBoxes = repmat([xCenter yCenter], data.stimuli.nSamples, 2) ...
    + [(posX-gaborDim) (posY-gaborDim) (posX+gaborDim) (posY+gaborDim)];

%----------------
% Run Experiment
%----------------

% Display instructions
DrawFormattedText(window, ['Welcome to the experiment! You will see an array of gratings.\n'...
    'Press the LEFT KEY if more gratings are orientated Clockwise.\n'...
    'Press the RIGHT KEY if more gratings are oriented Counter-Clockwise.\n'...,
    'Press any key to continue...'],...
    'center','center', white);

Screen('Flip', window);
KbStrokeWait;


% Run practice trials
DrawFormattedText(window, [sprintf('You will now see %i practice trials.\n', data.exp.pTrials)...
    'Press any key to continue...'],...
    'center','center', white);

Screen('Flip', window);
KbStrokeWait;

for i=1:data.exp.pTrials
    data = RunTrial(data, 0); % 0 -> no save data  
end

% Run blocks
for block=1:data.exp.numBlocks   
    for trial=1:data.exp.numTrialsPerBlock
        data = RunTrial(data, data.exp.numTrialsPerBlock*(block-1) + trial); % ~0 -> save data   
    end % end trials
   
    % Show block end info
    DrawFormattedText(window, [sprintf('You have finished %i blocks\n', block)...
        'You may choose to take a short break.\nPress any key to continue...'],...
        'center','center', white);
    
    Screen('Flip', window);
    
    KbStrokeWait;
    
end % end blocks

"trials ran"

% Save Data
SaveExit()

    % Function: Save and Exit
    function [] = SaveExit()
        Screen('CloseAll');
        sca;
        Priority(0);
        %save(fileName, 'data') 
        success = 1;
    end
end


%-----------------------------
% Functions from Target Paper
%-----------------------------

function [patch] = CreateCircularApertureINT(siz,falloff)
    if nargin < 2, falloff = 2; end
    if nargin < 1, error('Not enough input arguments.'); end
    
    sigmoidfun = @(x,lims)lims(1)+diff(lims)./(1+exp(-x));
    
    [x,y] = meshgrid((1:siz)-round((siz)/2));
    
    coef = log(1/0.01-1)*2/falloff;
    
    patch = sigmoidfun(coef*(sqrt(x.^2+y.^2)-siz/2),[1,0]);
end

