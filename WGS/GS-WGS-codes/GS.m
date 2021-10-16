clear all; close all;
tic;
%% Generate the SLM plane intensity: our resolution is 1920*1200
	x = linspace(-10,10,1920);
	y = linspace(-10,10,1200);
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
    Blank_pic = zeros(1200,1920);
    
%     Intensity_mat = [[255 255 255 255 255 255 255 255 255 255];
%                      [255 255 255 255 255 255 255 255 255 255];
%                      [255 255 255 255 255 255 255 255 255 255];
%                      [255 255 255 255 255 255 255 255 255 255];
%                      [255 255 255 255 255 255 255 255 255 255];
%                      [255 255 255 255 255 255 255 255 255 255];
%                      [255 255 255 255 255 255 255 255 255 255];
%                      [255 255 255 255 255 255 255 255 255 255];
%                      [255 255 255 255 255 255 255 255 255 255];
%                      [255 255 255 255 255 255 255 255 255 255]];
for i = 1:10
    for j = 1:10
        Blank_pic(800+20*i,700-30*j) = 80 + 4*abs(i-5)^2 + 4*abs(j-5.5)^2;%Intensity_mat(i,j);
    end
end

%% Prepare the target picture
	Target=double(Blank_pic);
	A = fftshift(ifft2(fftshift(Target)));
	error = [];
	iteration_num = 150;

%% Start Gershberg-Saxton algorithm iteration
for i=1:iteration_num
  B = abs(input_intensity) .* exp(1i*angle(A));
  C = fftshift(fft2(fftshift(B)));
  D = abs(Target) .* exp(1i*angle(C));
  A = fftshift(ifft2(fftshift(D)));
  error = [error; sum(sum(abs(1.32*abs(C) - abs(Target))))];   
end

%% Plot results
% 	figure
% 	subplot(2,1,1);
% 	imshow(Target);
% 	title('Original image')
% 	subplot(2,1,2);
% 	imagesc(abs(C))               % last pattern
% 	title('reconstructed image');
% 	figure
% 	i = 1:1:i;
% 	plot(i,(error'));
% 	title('Error');
% 	figure
% 	imagesc(abs(C)) %last pattern
% 	title('reconstructed image');
    SLM_phase = angle(A);
	imwrite(SLM_phase,'phase.png');
    imwrite(Target,"Target.png")
	toc;