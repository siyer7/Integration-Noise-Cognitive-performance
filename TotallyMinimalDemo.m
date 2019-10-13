% Clear the workspace
try
    close all;
    clearvars;
    sca;

    % Setup PTB with some default values
    PsychDefaultSetup(2);

    % Seed the random number generator. Here we use the an older way to be
    % compatible with older systems. Newer syntax would be rng('shuffle'). Look
    % at the help function of rand "help rand" for more information
    rand('seed', sum(100 * clock));

    % Set the screen number to the external secondary monitor if there is one
    % connected
    screenNumber = max(Screen('Screens'));

    % Define black, white and grey
    white = WhiteIndex(screenNumber);
    grey = white / 2;
    black = BlackIndex(screenNumber);

    % Open the screen
    [window, windowRect] = PsychImaging('OpenWindow', screenNumber, grey, [], 32, 2,...
        [], [],  kPsychNeed32BPCFloat);

    % Flip to clear
    Screen('Flip', window);

    % Query the frame duration
    ifi = Screen('GetFlipInterval', window);

    % Set the text size
    Screen('TextSize', window, 40);

    % Query the maximum priority level
    topPriorityLevel = MaxPriority(window);

    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(windowRect);
    w=xCenter; h=yCenter;
    [screenXpixels, screenYpixels] = Screen('WindowSize', window);
   %--------------------
% Gabor information
%--------------------

    % Dimension of the region where will draw the Gabor in pixels
    gaborDimPix = 100;%windowRect(4)/2;

    % Sigma of Gaussian
    sigma = gaborDimPix / 5;

    % Obvious Parameters
    orientation = 0;
    contrast = .6;
    aspectRatio = 1.0;
    phase = 0;

    % Spatial Frequency (Cycles Per Pixel)
    % One Cycle = Grey-Black-Grey-White-Grey i.e. One Black and One White Lobe
    numCycles = 5;
    freq = numCycles / gaborDimPix;

    % Build a procedural gabor texture (Note: to get a "standard" Gabor patch
    % we set a grey background offset, disable normalisation, and set a
    % pre-contrast multiplier of 0.5.
    % For full details see:
    % https://groups.yahoo.com/neo/groups/psychtoolbox/conversations/topics/9174
    backgroundOffset = [0.5 0.5 0.5 0.0];
    disableNorm = 1;
    preContrastMultiplier = 0.5;
    gabortex = CreateProceduralGabor(window, gaborDimPix, gaborDimPix, [],...
        backgroundOffset, disableNorm, preContrastMultiplier);

    % Randomise the phase of the Gabors and make a properties matrix.
    propertiesMat = [phase, freq, sigma, contrast, aspectRatio, 0, 0, 0];

    %----------------------------------------------------------------------
    %                       Timing Information
    %----------------------------------------------------------------------

    
    % fixation point time in seconds and frames
    fixationSecs = .3;
    fixationFrames = round(fixationSecs / ifi);

    % cue time in seconds and frames
    cueSecs = .7;
    cueFrames = round(cueSecs / ifi);
    
    % Stimulus array Time in seconds and frames
    gaborSecs = 0.15;
    gaborFrames = round(gaborSecs / ifi);

    % response time in seconds and frames
    responseSecs = 3;
    responseFrames = round(responseSecs / ifi);
    
    % feedback time in seconds and frames
    feedbackSecs = .5;
    feedbackFrames = round(feedbackSecs / ifi);
    
    % Numer of frames to wait before re-drawing
    waitframes = 1;


    %----------------------------------------------------------------------
    %                       Keyboard information
    %----------------------------------------------------------------------

    % Define the keyboard keys that are listened for. We will be using the left
    % and right arrow keys as response keys for the task and the escape key as
    % a exit/reset key
    escapeKey = KbName('ESCAPE');
    leftKey = KbName('LeftArrow');
    rightKey = KbName('RightArrow');


    %----------------------------------------------------------------------
    %                       Experimental loop
    %----------------------------------------------------------------------

    % Animation loop: we loop for the total number of trials
    r=400;angle=0;
    
    orientarr= [1 2 3 4 5 6 7 8];
    if(rand>0.5)
            avg= normrnd(3,8);
        else
            avg= normrnd(-3,8);
    end
    "avg "+avg
    for i = 1:8
        orientarr(i)= normrnd(avg,4);
        i+"= "+orientarr(i)
    end

    for trial = 1:5
        % Change the blend function to draw an antialiased fixation point
        % in the centre of the screen
        Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

        % If this is the first trial we present a start screen and wait for a
        % key-press
        if trial == 1
            DrawFormattedText(window, 'Press Any Key To Begin', 'center', 'center', black);
            Screen('Flip', window);
            KbStrokeWait;
        end

        % Flip again to sync us to the vertical retrace at the same time as
        % drawing our fixation point
        for frame=1:fixationFrames
            Screen('DrawDots', window, [xCenter; yCenter], 10, black, [], 2);
            vbl = Screen('Flip', window);
        end
        
        % cue
        for frame=1:cueFrames
            DrawFormattedText(window, 'Cue!', 'center', 'center', black);
            vbl = Screen('Flip', window);
        end
        
        % Gabor patches
        for frame = 1:gaborFrames
        
            DrawFormattedText(window, 'Cue!', 'center', 'center', black);
        
            % Set the right blend function for drawing the gabors
            Screen('BlendFunction', window, 'GL_ONE', 'GL_ZERO');

            for j=1:8
                Screen('DrawTextures', window, gabortex, [], [w+r*sin(angle)-50,h+r*cos(angle)-50,w+r*sin(angle)+50,h+r*cos(angle)+50],orientarr(j), [], [], [], [],...
                        kPsychDontDoRotation, propertiesMat');
                angle= angle+pi/4;
            end
            % Change the blend function to draw an antialiased fixation point
            % in the centre of the array
            Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');
            
            % Flip to the screen
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        end
        
        
        % response
%        for frame=1:responseFrames
        frame=1; respToBeMade = true;
        while(frame<=responseFrames && respToBeMade)
            frame=frame+1;
%            responseFrames
            DrawFormattedText(window, 'Cue!', 'center', 'center', black);
            DrawFormattedText(window, 'Response!', 'center', 100, black);
            vbl = Screen('Flip', window);
%            while respToBeMade
                [keyIsDown,secs, keyCode] = KbCheck;
                if keyCode(escapeKey)
                    ShowCursor;
                    sca;
                    return
                elseif keyCode(leftKey)
                    response = 1;
                    respToBeMade = false;
%                    DrawFormattedText(window, 'Left!', 'center', 'center'+200, black);
                elseif keyCode(rightKey)
                    response = 0;
                    respToBeMade = false;
%                    DrawFormattedText(window, 'Right!', 'center', 'center'+200, black);
                end
%               vbl = Screen('Flip', window);

%            end
        end
        
        % Now we present the isi interval with fixation point minus one frame
        % because we presented the fixation point once already when getting a
        % time stamp
        for frame = 1:feedbackFrames - 1

            % Draw the fixation point
            %Screen('DrawDots', window, [xCenter; yCenter], 10, black, [], 2);
            if respToBeMade
                DrawFormattedText(window, 'LATE!', 'center', 'center', black);
            else
                DrawFormattedText(window, 'Your Feedback!', 'center', 'center', black);
            end
            % Flip to the screen
            vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);
        end

        % Change the blend function to draw an antialiased fixation point
        Screen('BlendFunction', window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');

        % Flip to the screen
        vbl = Screen('Flip', window, vbl + (waitframes - 0.5) * ifi);

        % Now we wait for a keyboard button signaling the observers response.
        % The left arrow key signals a "left" response and the right arrow key
        % a "right" response. You can also press escape if you want to exit the
        % program
        
        
        % Record the response
        %respVector(stimValues == theAngle) = respVector(stimValues == theAngle)...
         %   + response;

        % Add one to the counter for that stimulus
        %countVector(stimValues == theAngle) = countVector(stimValues == theAngle) + 1;

    end

    %data = [stimValues; respVector; countVector]';

    %figure;
    %plot(data(:, 1), data(:, 2) ./ data(:, 3), 'ro-', 'MarkerFaceColor', 'r');
    %axis([min(data(:, 1)) max(data(:, 1)) 0 1]);
    %xlabel('Angle of Orientation (Degrees)');
    %ylabel('Performance');
    %title('Psychometric function');

    % Clean up
    sca;
catch
    Screen('CloseAll');
    %rethrow( ERR ); %now that the window is closed
end