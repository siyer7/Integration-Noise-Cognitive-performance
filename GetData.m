function [data] = GetData(varargin)
% Get Subject Data, load experimental parameters

data = struct;
data.participant = GetParticipantDataINT();

% Experimental Paramters
data.exp.numBlocks = 1;
data.exp.numTrialsPerBlock = 10;
data.exp.numTrials = data.exp.numBlocks * data.exp.numTrialsPerBlock;

data.exp.stimShowTime = 0.15; % 150ms
data.exp.responseTimout = 3; % 3sec, response timeout

% Keyboard Bindings
data.exp.escapeKey = KbName('ESCAPE');
data.exp.leftKey = KbName('LeftArrow');
data.exp.rightKey = KbName('RightArrow');

% Stimuli Parameters
data.stimuli.genMean = [-3 3];
data.stimuli.genStd = 8;
  
data.stimuli.contrast = [.15 .6];
data.stimuli.variability = [0 4 10];
data.stimuli.aspectRatio = 1.0;
  
data.stimuli.nSamples = 8; % number of gabors
data.stimuli.gaborDimPix = 100;    
data.stimuli.arrDiam = 200; % center screen -> center gabors

% Response Arrays
data.response.responseRight = Nan(1, data.exp.numTrials);
data.response.correct = Nan(1, data.exp.numTrials);
data.response.accuracy = Nan(1, data.exp.numTrials);
data.response.reactionTime = Nan(1, data.exp.numTrials);
data.response.cue = Nan(1, data.exp.numTrials);
data.response.orientationMean = Nan(1, data.exp.numTrials);
data.response.contrast = Nan(1, data.exp.numTrials);
data.response.variance = Nan(1, data.exp.numTrials);
data.response.trueOrientaions = Nan(data.stimuli.nSamples, data.exp.numTrials);



end

function [participant] = GetParticipantDataINT()

    argindlg = inputdlg({'Participant ID','Gender (M/F/X)','Age','Hand (L/R)'}, ...
        '',1,{'SUB','','','R'});
    if isempty(argindlg)
        participant = struct;
        participant.name = 'NULL';
    else
        participant = struct;
        participant.name = ['B' upper(argindlg{1})];
        participant.gender = argindlg{2};
        participant.age  = argindlg{3};
        participant.handedness = upper(argindlg{4});
        participant.now = datestr(now,'yyyymmddHHMMSS');
    end

end