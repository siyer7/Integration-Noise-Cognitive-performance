Here we provide scripts to the audio recordings in each section of our poster.
Recordings indicated by section title.

## Introduction
Are you confident in your decisions? Do you think you are justified in your confidence or could you potentially be overconfident? Previous findings do not quite agree on whether we are justified in our confidence. For our study, we started with a replication which found that when people are faced with a certain kind of noise their confidence tends to be less aligned with their accuracy. We looked at three different conditions which manipulated the type of noise displayed to subjects: one condition which had no noise, one which had encoding noise (making it more difficult to see), and one which had integration noise (making it more difficult to process). We predicted that people’s confidence would be misaligned when faced with integration noise. However, our results were mixed, which warrants future work. For more information about our study, feel free to click on the audio or script for each section on our poster.

### Intuition about Encoding and Integration Noise
To build an intuition of these noises, let me walk you through a short example. Imagine that your friend is throwing a ball at you. If it is dark outside, you might be able to encode this information, and confidently and accurately adjust your position to catch the ball. On the other hand, if it is windy, you might not integrate this information as well, and might confidently but inaccurately adjust your position, since the ball might land far from you. Encoding noise refers to it being dark outside; integration noise refers to it being windy outside, in which case you are confident but inaccurate in your position. Note that this example is mainly for illustration purposes, and that these terminologies and definitions are not agreed upon in the field of cognitive science.

## Methods
Our experiment had 3 conditions. We increase the integration noise by increasing the variability of orientations across the patches. We increase the encoding noise by decreasing the contrast of the patches. To recap, the high variability condition presents high integration noise, while the low contrast condition presents high encoding noise. A 3rd condition, called the baseline, serves as a control. To measure subjects’ confidence, we study their reliance on a hint that we provide, and also ask them to explicitly report their confidence. Presumably, if a subject has high high confidence in their performance, or accuracy, they will show lesser reliance on our cue, and express higher self-reported confidence.

# Results

## Lower accuracy in 'high variability'
We want to know whether our stimulus manipulation did in fact add noise, 
that is make the task more difficult. We compared subjects’ accuracy across
the different conditions and found that, on average, performance was highest
in the ‘baseline’ and lowest in ‘high variance’. However, there are 
individual differences in this trend. Some subjects (like subject 4), 
performed almost as well in ‘low contrast’ as in ‘baseline’. The inline 
plot shows the difference in accuracy between ‘low contrast’ and ‘high 
variability’. Notice that almost all subjects performed better in ‘low 
contrast’ than in ‘high variability’.

## Lower Cue Reliance in 'high variability'
In the experiment, we provided an informative cue on half the trials. We 
assume subjects will rely more on the cue if they are less confident, so our
hypothesis predicts that subjects will be less influenced by the cue in the
‘high variability condition’.
To calculate how much subjects relied on the cue, we fit psychometric curves
to subjects’ responses for each of the cues presented in each condition. 
Psychometric curves model how a subject’s response changes as a function of
the signal in the stimulus (in our case, as the average orientation moves
away from 0, it becomes easier to determine the correct answer). Looking at
the example provided, we see that this subject’s responses were indeed 
influenced by the cue, as seen by the shift in curves in the presence of 
cues. To quantify this, we calculated the distance between these curves.
We found that, on average, subjects were indeed less influenced by the cue
in the ‘high variance condition’. Individual subjects showed some differences.
Looking at the inline plot which shows the difference in reliance between
the ‘low contrast’ and ‘high variance’ conditions, we see that no subjects
were more reliant on cues in the ‘high variability’ condition, giving us 
weak evidence that people are less able to account for integration noise
in their confidence.

## Similar confidence-accuracy alignment in ‘high variability’ and ‘low contrast’
Besides cue usage, we also analyzed subjects’ self-reported confidence levels. We hypothesized a greater misalignment between confidence & accuracy in the high variability, that is, high integration noise, condition. To this end, we compute the alignment, or spearman correlation, between confidence and accuracy. We then bootstrap - in simple terms, this means that we repeat this procedure a 1000 times, to create a probability distribution and sharpen our evidence. Next, we compare this correlation between the low contrast & high variability conditions. In the plot, dotted & color lines denote subject behavior, while the black line denotes the combined population behavior. We see that the difference between the 2 conditions, denoted by the x-axis, is centered very close to 0, denoted by the vertical line. This suggests that there is no significant difference between the 2 conditions. Had the black plot been centered at a positive number, instead of 0, it would suggest a misalignment in the high variability condition. In the smaller bar plot, we see that some subjects had their values centered at a positive number, and some at a negative number. Thus, we have a mixed bag of results, which does not support or reject our hypothesis with clarity.

## Discussion
Overall, we found that subjects, on average, rely about the same or less
on cues in the ‘high variability’ condition, despite having similar or
worse accuracy, giving us weak evidence that people cannot account for
integration noise in their confidence. However, analysis of the explicit
confidence measure found mixed results. On average, subjects did not have
better aligned confidence in the ‘low contrast’ condition. Thus we did 
not successfully replicate the target paper.
Extending this work by looking at individual subjects, we see more mixed
effects. Some subjects were clearly less reliant on cues in the high 
variance condition, and some subjects clearly had better aligned confidence
in ‘low contrast’, as predicted by our hypothesis. However, some subjects
show an opposite trend, wherein they had better aligned confidence in the
‘high variance’ condition, the opposite of what we predict.
One of the biggest limitations of our study was that subjects’ accuracies
were not the same in the two experimental conditions. We expect that, in
general, subjects will be more confident if they are performing better.
This means that we cannot tease apart smaller confidence from worse 
performance. To resolve this, we would repeat the experiment using a
staircase method, which will ensure that each subject is performing with
the same accuracy in the experimental conditions. Furthermore, subjects
testing was cut short due to COVID-19, so the statistical power of our
experiment is lower than we anticipated.



