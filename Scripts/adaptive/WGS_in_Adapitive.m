clear all; close all;
tic;
%% Parameters

row = 10;
column = 10;
row_spacing = 15;
column_spacing = 10;
iteration_num = 20;
Resolution = [1920 1200];
adaptive_num = 50;
Camera_Exposure = 300;%112.181;
Camera_Gain = 1;


%% Initialize the camera
% create a video obj to get data
vidobj = videoinput('gentl', 1, 'Mono8'); 
vidobj.TriggerRepeat=inf;
vidobj.FramesPerTrigger=1;
% set trigger mode to manual, then we can start and stop the
% camera by our self
triggerconfig(vidobj, 'manual');
% Get the video source
src = getselectedsource(vidobj); 
start(vidobj);

%% Generate the SLM plane intensity: our resolution is Resolution(1)*Resolution(2)
	x = linspace(-10,10,Resolution(1));
	y = linspace(-10,10,Resolution(2));
	[X,Y] = meshgrid(x,y);
	x0 = 0;     		% center
	y0 = 0;     		% center
	sigma = 5; 			% effective beam waist
    sigma_x = sigma/1.92;
    sigma_y = sigma/1.2;
	A = 1;      		% peak of the beam 
	res = ((X-x0).^2./(2*sigma_x^2) + (Y-y0).^2./(2*sigma_y^2));
	input_intensity = A  * exp(-res);
% 	surf(input_intensity);
% 	shading interp


%% Generate the first target picture
    Blank_pic = zeros(Resolution(2),Resolution(1));

    for i = 1:column
        for j = 1:row
            Blank_pic(600+column_spacing*i,960-row_spacing*j) = 1.0;
        end
    end

    % Prepare the target picture
	Target=double(Blank_pic);
	A = fftshift(ifft2(fftshift(Target)));

%% LOOP starts
% record = [];
stdeviation = [];

for adaptive = 1:adaptive_num
%% Start Weighted Gershberg-Saxton algorithm iteration
	weight = ones(column,row);
    for k=1:iteration_num
        B = abs(input_intensity) .* exp(1i*angle(A));
        B = B / max(max(abs(B)));
        C = fftshift(fft2(fftshift(B)));
        amplitute_retrival = abs(C);
        aver_intensity_retrival = 0;
    
        for i = 1:column % calculate retrival intensity, next we can weighted
            for j = 1:row
                aver_intensity_retrival = amplitute_retrival(600+column_spacing*i,960-row_spacing*j) + aver_intensity_retrival;
            end
        end
        aver_intensity_retrival = aver_intensity_retrival / row / column;
        
        for i = 1:column % weighted part
            for j = 1:row
                weight(i,j) = weight(i,j) / (amplitute_retrival(600+column_spacing*i,960-row_spacing*j) / aver_intensity_retrival);
            end
        end
    
        D = abs(Target) .* exp(1i*angle(C));
        
        for i = 1:column % weighted part
            for j = 1:row
                D(600+column_spacing*i,960-row_spacing*j) = D(600+column_spacing*i,960-row_spacing*j) * weight(i,j);
            end
        end
        
        D = D/max(max(abs(D)));
        A = fftshift(ifft2(fftshift(D)));   
    end
    pause(0.1);

    SLM_phase = im2uint8((angle(A)+pi)/2/pi);
% 	imwrite(SLM_phase,'phase.png');
    pause(0.1);
%% Display into second monitor

    Screen('Preference', 'SkipSyncTests', 2);
    window = Screen('OpenWindow',2);
    Screen('PutImage',window,SLM_phase);
    Screen('Flip',window); 

    pause(0.1);
%% Take photo

    src.ExposureTime = Camera_Exposure;
    src.Gain = Camera_Gain;
    snapshot = getsnapshot(vidobj);
    pause(0.1);
    snapshot = getsnapshot(vidobj);
    pause(0.1);
    
%% Analyse tweezer info

    BW = imbinarize(snapshot);
    [B,~] = bwboundaries(BW,'noholes');

    tweezer_information = zeros(length(B),1); %Store tweezer position and intensity
    for k = 1:length(B) % for k-th tweezer
       boundary = B{k};
       max_x = max(boundary(:,2));
       min_x = min(boundary(:,2));
       min_y = min(boundary(:,1));
       max_y = max(boundary(:,1));

       current_tweezer = snapshot(min_y:max_y,min_x:max_x);
%        tweezer_information(k) = sum(sum(current_tweezer));
       max_value = max(max(current_tweezer));
       tweezer_information(k) = max_value;
    end

    % the tweeze information is currently intensity
    tweezer_number = row * column;
    if tweezer_number ~= length(tweezer_information)
        tweezer_information(1) = []; % get rid of 0-order
    end

    stdeviation = [stdeviation; std(tweezer_information)/mean(tweezer_information)];
    tweezer_information = reshape(tweezer_information,column,row);
    tweezer_information = double(tweezer_information);
    tweezer_information = tweezer_information / mean(mean(tweezer_information));
    
%     record = [record; (max(max(tweezer_information))-min(min(tweezer_information)))/(max(max(tweezer_information))+min(min(tweezer_information)))];
%% Update the target figure

    for i = 1:column
        for j = 1:row
            Blank_pic(600+column_spacing*i,960-row_spacing*j) = Blank_pic(600+column_spacing*i,960-row_spacing*j) / sqrt(tweezer_information(i,j));
        end
    end
    % Prepare the target picture
	Target=double(Blank_pic).* exp(1i*angle(C));
	A = fftshift(ifft2(fftshift(Target)));

% Show figure
adaptive_list = 1:1:adaptive;
% plot(adaptive_list, record);
plot(adaptive_list, stdeviation);
end


%% Shut down the camera
    stop(vidobj);
    delete(vidobj);
%% Show the last figure
	toc;