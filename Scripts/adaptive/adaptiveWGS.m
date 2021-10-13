clear all; close all;
tic;
%% Parameters

row = 5;
column = 10;
row_spacing = 15;
column_spacing = 10;
iteration_num = 30;
Resolution = [1920 1200];
adaptive_num = 5;
Camera_Exposure = 112.181;
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

    % % set the Binning of the camera
    % src.BinningHorizontal = 2;
    % src.BinningVertical = 2;

    % start video aquisition
    start(vidobj);
    
    % get a pic
    trigger(vidobj);
    
    %snapshot = getsnapshot(vidobj); 
    snapshot=getdata(vidobj);
    
    % adjust to the size of our pic
%     im = image(zeros(size(snapshot),'uint8'));
    
    % set Preview as image mode
    axis('image');

    
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

    weight = ones(column,row);

    % Prepare the target picture
	Target=double(Blank_pic);
	A = fftshift(ifft2(fftshift(Target)));


    

%% LOOP starts
record = [];

for adaptive = 1:adaptive_num
%% Start Weighted Gershberg-Saxton algorithm iteration
	error = [];
    for k=1:iteration_num
        B = abs(input_intensity) .* exp(1i*angle(A));
        B = B / max(max(abs(B)));
        C = fftshift(fft2(fftshift(B)));
        intensity_retrival = abs(C);
        aver_intensity_retrival = 0;
    
        for i = 1:column % calculate retrival intensity, next we can weighted
            for j = 1:row
                aver_intensity_retrival = intensity_retrival(600+column_spacing*i,960-row_spacing*j) + aver_intensity_retrival;
            end
        end
        aver_intensity_retrival = aver_intensity_retrival / row / column;
        
        for i = 1:column % weighted part
            for j = 1:row
                weight(i,j) = weight(i,j) / (intensity_retrival(600+column_spacing*i,960-row_spacing*j)/aver_intensity_retrival);
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
        error = [error; sum(sum(abs(intensity_retrival/sum(sum(intensity_retrival)) - abs(Blank_pic))))];   
    end
    pause(0.1);
%% Plot results
% 	figure
% 	subplot(2,1,1);
% 	imshow(Target);
% 	title('Original image')
% 	subplot(2,1,2);
% 	imagesc(abs(C))               % last pattern
% 	title('reconstructed image');
% 	figure
% 	k = 1:1:k;
% 	plot(k,(error'));
% 	title('Error');
% 	figure
% 	imagesc(abs(C)) %last pattern
% 	title('reconstructed image');
    SLM_phase = im2uint8((angle(A)+pi)/2/pi);
	imwrite(SLM_phase,'phase.png');
%     imwrite(Target,"Target.png")
%     imwrite(abs(C),"Retrived.png")
    pause(0.1);
%% Display into second monitor

    Screen('Preference', 'SkipSyncTests', 2);
    window = Screen('OpenWindow',2);
%     image = imread("C:\Users\doyle\Documents\SLM-Tweezer\3by3_tweezer\phase.png");
%     Screen('PutImage',window,image);
    Screen('PutImage',window,SLM_phase);
    Screen('Flip',window); 

    pause(0.1);
%% Take photo

    src.ExposureTime = Camera_Exposure;
    src.Gain = Camera_Gain;
    snapshot = getsnapshot(vidobj);
    pause(0.3);
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
       max_value = max(max(current_tweezer));
    %    [max_value_y,max_value_x] = find(max_value == current_tweezer);
    %    max_value_y = max_value_y(1) + min_y - 1;
    %    max_value_x = max_value_x(1) + min_x - 1;
    %    tweezer_information{k} = {max_value_y; max_value_x; max_value};
       tweezer_information(k) = max_value;
    end

    % the tweeze information is currently position y and x and value
    tweezer_number = row * column;
    if tweezer_number ~= length(tweezer_information)
        tweezer_information(1) = [];
    end
    tweezer_information = reshape(tweezer_information,column,row);
    tweezer_information = double(tweezer_information);
    tweezer_information = tweezer_information / mean(mean(tweezer_information));
    
    record = [record; (max(max(tweezer_information))-min(min(tweezer_information)))/(max(max(tweezer_information))+min(min(tweezer_information)))];
%% Update the target figure

    for i = 1:column
        for j = 1:row
            Blank_pic(600+column_spacing*i,960-row_spacing*j) = Blank_pic(600+column_spacing*i,960-row_spacing*j) / tweezer_information(i,j);
        end
    end

    weight = ones(column,row);

    % Prepare the target picture
	Target=double(Blank_pic);
	A = fftshift(ifft2(fftshift(Target)));

%% Another way to update -- not right!
%     for i = 1:column
%         for j = 1:row
%           weight(column,row) = weight(column,row) / tweezer_information(i,j);
%         end
%     end
%% Loop ends
adaptive_list = 1:1:adaptive;
plot(adaptive_list, record);
end


%% Shut down the camera
    % we must stop this video object by ourself, otherwise problem
    % will come
    stop(vidobj);

    % delete the video object, so other thread can use this camera
    delete(vidobj);
%% Show the last figure
	toc;