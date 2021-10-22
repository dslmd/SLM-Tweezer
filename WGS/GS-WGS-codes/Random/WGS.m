clear all; close all;
tic;

row = 10;
column = 10;
row_spacing = 15;
column_spacing = 10;
iteration_num = 30;
Resolution = [1920 1200];
tweezer_number = row * column;

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

%% Generate the target picture
    Blank_pic = zeros(Resolution(2),Resolution(1));

    for i = 1:column
        for j = 1:row
            Blank_pic(600+column_spacing*i,960-row_spacing*j) = 1;
        end
    end

    weight = ones(column,row);

    % Prepare the target picture
	Target=double(Blank_pic);
	A = fftshift(ifft2(fftshift(Target)));
	error = [];

%% Start Weighted Gershberg-Saxton algorithm iteration
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
%         error = [error; sum(sum(abs(intensity_retrival/sum(sum(intensity_retrival)) - abs(Blank_pic))))];   
    end
    non_uniformity = [];
% for sigma = 0.1:0.1:5
% %% Intensity noise?
% 	x = linspace(-10,10,Resolution(1));
% 	y = linspace(-10,10,Resolution(2));
% 	[X,Y] = meshgrid(x,y);
% 	x0 = 5;     		% center
% 	y0 = 0;     		% center
%  	sigma = 5; 			% effective beam waist
%     sigma_x = sigma/1.92;
%     sigma_y = sigma/1.2;
% 	res = ((X-x0).^2./(2*sigma_x^2) + (Y-y0).^2./(2*sigma_y^2));
% 	input_intensity = exp(-res);

%% final loop, add noise to see retrive
    % add phase noise in B
    B = abs(input_intensity) .* exp(1i*angle(A)) .* exp(1i * ((rand(1200,1920)*2-1) * 0.1 * pi));
    B = B / max(max(abs(B)));
    C = fftshift(fft2(fftshift(B)));

%% Plot results

    Retrive = abs(C)/max(max(abs(C)));
    SLM_phase = im2uint8((angle(A)+pi)/2/pi);
	imwrite(SLM_phase,'phase.png');
    imwrite(Target,"Target.png")
    imwrite(Retrive,"Retrived.png")
	
    
    intensity_distribution = zeros(column,row);
    for i = 1:column % weighted part
        for j = 1:row
            intensity_distribution(i,j) = Retrive(600+column_spacing*i,960-row_spacing*j);
        end
    end

    non_uniformity = [non_uniformity; std(reshape(intensity_distribution,100,1))/mean(mean(intensity_distribution))]

% end
% plot(non_uniformity);
toc;
    