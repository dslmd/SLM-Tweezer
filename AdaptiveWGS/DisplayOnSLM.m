% We first need to initialize the toolbox.
Screen('Preference', 'SkipSyncTests', 2);
image = imread("C:\Users\doyle\Documents\SLM-Tweezer\3by3_tweezer\phase.png");
window = Screen('OpenWindow',2);
Screen('PutImage',window,image);
Screen('Flip',window); 