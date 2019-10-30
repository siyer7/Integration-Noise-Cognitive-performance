function [data] = Experiment(varargin)
% Main Code

% Clear Workspace
close all;
sca;
clear vars;

%% File handling
directory = fullfile(pwd, 'data');

% check existence (and make) directory
if ~exist(directory, 'dir')
    mkdir(directory);
end

% Get Subject Data and Parameters
data = GetData(directory);
 
fileName = fullfile(directory, [data.participant.subjectId '-1.mat']);
fileNameEarly = fullfile(directory, [data.participant.subjectId '-1-EarlyQuit' '.mat']);

i = 1;
while exist(fileName, 'file') || exist(fileNameEarly, 'file')
    i = i+1;
    fileName = fullfile(directory, [data.participant.subjectId sprintf('-%i', i) '.mat']);
    fileNameEarly = fullfile(directory, [data.participant.subjectId sprintf('-%i', i) '-EarlyQuit' '.mat']);
end


%% Prepare Psychtoolbox Stuff
PsychDefaultSetup(2);

screenNum = max(Screen('Screens'));

black = [ 0  0  0];
grey =  [.5 .5 .5];
white = [ 1  1  1];

[window, windowRect] = PsychImaging('OpenWindow', screenNum, grey, [1920 0 1920*2 1080], 32, 2,...
        [], [],  kPsychNeed32BPCFloat);

[xCenter, yCenter] = RectCenter(windowRect);
[xWidth, yHeight] = Screen('WindowSize', window); 
Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

ifi = Screen('GetFlipInterval', window);
stimShowTime = data.exp.stimShowTime;
data.exp.stimFrames = round(stimShowTime/ifi);

topPriorityLevel = MaxPriority(window);
Priority(topPriorityLevel);
HideCursor;

Screen('TextSize', window, data.exp.directionTextSize);

% Prepare Stimulus Stuff
data.stimuli.aperture = CreateCircularApertureINT(data.stimuli.gaborDimPix, 1);

data.stimuli.fixBox = [xCenter yCenter xCenter yCenter] + [-1 -1 1 1].*2; % location of fixation

posX = data.stimuli.posX'; 
posY = data.stimuli.posY';
gaborDim = data.stimuli.gaborDimPix;
data.stimuli.gratingBoxes = repmat([xCenter yCenter], data.stimuli.nSamples, 2) ...
    + [(posX-gaborDim) (posY-gaborDim) (posX+gaborDim) (posY+gaborDim)];


%% Run Experiment

% Display instructions
DrawFormattedText(window, ['Welcome to the experiment! You will see an array of gratings.\n'...
    'Press the 1 KEY if more gratings are orientated Clockwise.\n'...
    'Press the 3 KEY if more gratings are oriented Counter-Clockwise.\n'...,
    'Press any key to continue...'],...
    'center','center', white);

Screen('Flip', window);
KbStrokeWait;
  

% Inform about practice trials
DrawFormattedText(window, [sprintf('You will now see %i practice trials.\n', data.exp.pTrials)...
    'Press any key to continue...'],...
    'center','center', white);


Screen('Flip', window);
KbStrokeWait;

Screen('TextSize', window, data.exp.cueTextSize);

% Run practice trials
for i=1:data.exp.pTrials
    cueBlock= binornd(1,0.5); % 0.5 probability of block having or not having cues
    [data, timedOut, quit] = RunTrial(data, 0, window); % no save data 
    if quit
        Screen('CloseAll');
        sca;
        return
    end
end

% Pre-block text
Screen('TextSize', window, data.exp.directionTextSize);

DrawFormattedText(window, 'Are you ready for Block 1?\nPress any button to continue...',...
    'center','center', white);

Screen('Flip', window);
KbStrokeWait;

% Run blocks
block = 1;
count_uncued=0
count_cued=0
while block<=data.exp.numBlocks
    %%% ensuring no. of cued blocks= no.of uncued blocks
    prob= (18-count_cued)/36
    cueBlock= binornd(1,prob); % 0.5 probability of block having or not having cues
    if cueBlock==0
        count_uncued= count_uncued+1
    else
        count_cued= count_cued+1
    end
    trial = 1; % run trials
    while trial<=data.exp.numTrialsPerBlock
        data.response.isCuedBlock(:,data.exp.numTrialsPerBlock*(block-1) + trial) = cueBlock; % save whether block is cued
        [data, timedOut, quit] = RunTrial(data, data.exp.numTrialsPerBlock*(block-1) + trial, window); % save data
        
        if quit
            SaveExit(quit);
            return;
        end
        
        if ~timedOut
            trial = trial + 1;
        end
    end % end trials
    
    % Show block end info
    DrawFormattedText(window, [sprintf('You have finished %i blocks\n', block)...
        'You may choose to take a short break.\nPress any key to continue...'],...
        'center','center', white);
    
    Screen('Flip', window);
    
    KbStrokeWait;
    
    block = block + 1;
    
end % end blocks

"trials ran" % DEBUG

%% Save Data
SaveExit()

    % Save and Exit
    function [] = SaveExit(early)
        if nargin==0
            early=0;
        end
        
        Screen('CloseAll');
        sca;
        Priority(0);
        
        try % some trials have been performed
            numTrials = data.exp.numTrialsPerBlock*(block-1) + trial;
        catch % no trials --> don't save anything
            return
        end
        
        if early
            save(fileNameEarly, 'data');
        else
            save(fileName, 'data');
        end
        "Succesfully Saved!"
    end % end SaveExit

end % end Experiment


%% Functions from Target Paper

function [patch] = CreateCircularApertureINT(siz,falloff)
    if nargin < 2, falloff = 2; end
    if nargin < 1, error('Not enough input arguments.'); end
    
    sigmoidfun = @(x,lims)lims(1)+diff(lims)./(1+exp(-x));
    
    [x,y] = meshgrid((1:siz)-round((siz)/2));
    
    coef = log(1/0.01-1)*2/falloff;
    
    patch = sigmoidfun(coef*(sqrt(x.^2+y.^2)-siz/2),[1,0]);
end

