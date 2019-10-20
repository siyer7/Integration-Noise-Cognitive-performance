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


% Prepare Psychtoolbox
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
    data = RunTrial(data, 1); % 1 -> no save data  
end

% Run blocks
for block=1:data.exp.numBlocks   
    for trial=1:data.exp.numTrialsPerBlock
        data = RunTrial(data, 0); % 0 -> save data   
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


