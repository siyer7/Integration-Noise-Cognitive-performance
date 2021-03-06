function [data] = GetData(directory, varargin)
% Get Subject Data & Load Experimental Parameters
% Set all the experimental parameters here

data = struct;
data.participant.subjectId = getSubjectId(directory, 'uncertaintyV1')

% Experimental Paramters
data.exp.numBlocks = 25;
data.exp.numTrialsPerBlock = 36;
data.exp.numTrials = data.exp.numBlocks * data.exp.numTrialsPerBlock;

data.exp.pTrials = 72; % number of practice trials

data.exp.stimShowTime = .150; % 150ms
data.exp.responseTimout = 3; % 3sec, response timeout

data.exp.fixShowTime = 0.3; % 300ms show fixation point
data.exp.cueShowTime = 0.7; % 700ms show cue
data.exp.feedbackShowTime = 1; % 1sec show feedback

data.exp.rangSeeds = NaN(1, data.exp.numTrials);

data.exp.cueTextSize = 16;
data.exp.directionTextSize = 40;

% Keyboard Bindings
KbName('UnifyKeyNames');
data.exp.escapeKey = KbName('ESCAPE');
data.exp.leftKey = KbName('4');
data.exp.rightKey = KbName('6');
data.exp.high = KbName('8');
data.exp.medium= KbName('5');
data.exp.low = KbName('2');

% Stimuli Parameters
data.stimuli.genMean = [-3 3];
data.stimuli.genStd = 8;
  
data.stimuli.contrastVal = [.15 .6];
data.stimuli.variabilityVal = [0 10];
  
data.stimuli.nSamples = 8; % number of gabors
data.stimuli.gaborDimPix = 50;    
data.stimuli.arrDiam = 200; % center screen -> center gabors

data.stimuli.aspectRatio = 1.0;
data.stimuli.frequency = 0.05*2;
data.stimuli.envelopDev = 90;
data.stimuli.noiseAmp = .1;
data.stimuli.phase = pi*rand(data.stimuli.nSamples, data.exp.numTrials);

% Calculate Stimuli position for all 8 gabor patches
allAngles = linspace(0, 360, data.stimuli.nSamples+1); 
allAngles = pi ./ 180 .* allAngles(1:end-1);

posX = data.stimuli.arrDiam .* cos(allAngles); 
posY = data.stimuli.arrDiam .* sin(allAngles);

data.stimuli.posX = round(posX);
data.stimuli.posY = round(posY);

% Response Arrays
data.response.randSeed = NaN(1, data.exp.numTrials);

data.response.responseRight = NaN(1, data.exp.numTrials); % subject response
data.response.correct = NaN(1, data.exp.numTrials);
data.response.accuracy = NaN(1, data.exp.numTrials);
data.response.reactionTime = NaN(1, data.exp.numTrials);
data.response.confidence = NaN(1, data.exp.numTrials);

data.response.isCuedBlock = NaN(1, data.exp.numTrials); % boolean
data.response.cue = NaN(1, data.exp.numTrials);
data.response.orientationMean = NaN(1, data.exp.numTrials);
data.response.contrast = NaN(1, data.exp.numTrials);
data.response.variance = NaN(1, data.exp.numTrials);
data.response.trueOrientaions = NaN(data.stimuli.nSamples, data.exp.numTrials);

end

% Screen to get participant infos (NOT USED ANYMORE)
function [participant] = GetParticipantDataINT()

    argindlg = inputdlg({'Participant ID','Gender (M/F/X)','Age','Hand (L/R)'}, ...
        '',1,{'SUB','','','R'});
    if isempty(argindlg)
        participant = struct;
        participant.name = 'NULL';
    else
        participant = struct;
        participant.name = getSubjectId(argindlg{1});
        participant.gender = argindlg{2};
        participant.age  = argindlg{3};
        participant.handedness = upper(argindlg{4});
        participant.now = datestr(now,'yyyymmddHHMMSS');
    end

end