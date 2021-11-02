clear all;
SLM_Phase = ones(1200,1920);
mask_radius = 500;
C31 = -0.05;
C00 = 0.85 - abs(C31);


%% Mask the SLM phase
SLM_Phase = double(SLM_Phase);
SLM_Phase = SLM_Phase/max(max(SLM_Phase));
SLM_Phase = SLM_Phase*2*pi - pi;
pic_shape = size(SLM_Phase);
Masked_Phase = SLM_Phase;

%%
for i = 1:pic_shape(1)
    for j = 1:pic_shape(2)
        r = sqrt( (i-pic_shape(1)/2)^2 + (j-pic_shape(2)/2)^2 ) / mask_radius;
        if r <= 1
            phi = atan((j-pic_shape(2)/2)/(i-pic_shape(1)/2));
            Masked_Phase(i,j) = SLM_Phase(i,j) * (C00 * Zernike(0,0,r,phi) +  C31 * Zernike(2,2,r,phi) +  0.15 * Zernike(4,0,r,phi));
        end
    end
end
%%

Masked_Phase = (Masked_Phase + pi)/ 2 / pi;
imwrite(Masked_Phase,'Masked_Phase.png');