SLM_Phase = imread('/Users/dslmd/Downloads/SLM-Tweezer/WGS/Example/10by10_tweezer/Camera.bmp');
row = 10;
column = 10;
pixel_size = 1;

%% Start with an image

BW = imbinarize(I);
[B,~] = bwboundaries(BW,'noholes');
tweezer_information = cell(length(B),1); %Store tweezer position and intensity

for k = 1:length(B) % for k-th tweezer
    boundary = B{k};
    max_x = max(boundary(:,2));
    min_x = min(boundary(:,2));
    min_y = min(boundary(:,1));
    max_y = max(boundary(:,1));
    
    current_tweezer = I(min_y:max_y,min_x:max_x);

    % Here we are in the single tweezer region, so we can fit seperately
    % grab a new figure from camera
    max_value = max(max(current_tweezer));
    [y_max, x_max] = find(max_value == current_tweezer);
    
    % extract x data
    x = sum(current_tweezer,1);
    x_size = size(x);
    x_array = [1:1:x_size(2)];
    x_max_value = max(x);
    x = x';
    x_array = x_array';
    
    % extract y data
    y = sum(current_tweezer,2);
    y_size = size(y);
    y_array = [1:1:y_size];
    y_array = y_array';
    y_max_value = max(y);
    
    f = fittype('a*exp(-((x-b)/c)^2)+d','independent','x','coefficients',{'a','b','c','d'});
     
    % fit x
    c_fun_x = fit(x_array,x,f,'StartPoint',[x_max_value, x_max(1), 10,0],'Lower',[1,0,1,0],'Upper',[inf,inf,inf,inf]);
    x_max_value = c_fun_x.a * pixel_size;
    x_max = (min_x + c_fun_x.b - 1) * pixel_size;
    x_waist = c_fun_x.c * pixel_size;
    
    % fit y
    c_fun_y = fit(y_array,y,f,'StartPoint',[y_max_value, y_max(1), 10,0],'Lower',[1,0,1,0],'Upper',[inf,inf,inf,inf]);
    y_max_value = c_fun_y.a * pixel_size;
    y_max = (min_y + c_fun_y.b - 1) * pixel_size;
    y_waist = c_fun_y.c * pixel_size;

    tweezer_information{k} = [x_max_value, x_max, x_waist, y_max_value, y_max, y_waist];
end

%% the tweeze information is currently position y and x and value

tweezer_number = row * column;
if tweezer_number ~= length(tweezer_information)
    tweezer_information(1) = [];
end

% tweezer_information = reshape(tweezer_information,column,row);
% tweezer_information = tweezer_information';

intensity = zeros(1,row*column);
for i = 1:row*column
    intensity(i) = sqrt(tweezer_information{i}(1) * tweezer_information{i}(4));
end
non_uniformity = std(intensity)/mean(intensity);

astigmatism = zeros(1,row*column);
for i = 1:row*column
    astigmatism(i) = max(tweezer_information{i}(3)/tweezer_information{i}(6), tweezer_information{i}(6)/tweezer_information{i}(3));
end
max_astigmatism = max(astigmatism);