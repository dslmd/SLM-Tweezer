%% "I" is just our camera figure

I = imread('/Users/dslmd/Downloads/SLM-Tweezer/5by5_tweezer/Camera.bmp');
row = 5;
column = 5;

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
   max_value = max(max(current_tweezer));
%    [max_value_y,max_value_x] = find(max_value == current_tweezer);
%    max_value_y = max_value_y(1) + min_y - 1;
%    max_value_x = max_value_x(1) + min_x - 1;
%    tweezer_information{k} = {max_value_y; max_value_x; max_value};
   tweezer_information{k} = max_value;
end

%% the tweeze information is currently position y and x and value

tweezer_number = row * column;
if tweezer_number ~= length(tweezer_information)
    tweezer_information(1) = [];
end
tweezer_information = reshape(tweezer_information,column,row);
tweezer_information = tweezer_information';
