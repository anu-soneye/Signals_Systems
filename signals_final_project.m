
%% Examples With Windowing


f = [0 0.48 0.48 1];            % Frequency breakpoints
m = [0 0 1 1];                  % Magnitude breakpoints

b = fir1(34,0.48,'high',hamming(35));   % FIR filter design
freqz(b,1,512)

d = fir1(34,0.48,'high',rectwin(35));   % FIR filter design
[H3, w3] = freqz(d,1,512);          % Frequency response of filter
freqz(d,1,512);


%% Chirp Analysis

load chirp
close all
y2 = y + randn(size(y))/25;
plot(y2)
%sound(y3)
f = [0 0.3 0.3 1];            % Frequency breakpoints
m = [0 0 1 1];                  % Magnitude breakpoints
b = fir2(34,f,m,hamming(35));   % FIR filter design
c = fir2(200,f,m,hamming(201));   % FIR filter design

f2 = [0 0.65 0.65 1];            % Frequency breakpoints
m2 = [0 0 1 1];                  % Magnitude breakpoints
d = fir2(34,f2,m2,hamming(35));   % FIR filter design
e = fir2(34,f2,m2,chebwin(35, 27));   % FIR filter design

[H1, w1] = freqz(b,1,512);         % Frequency response of filter
%freqz(b,1,512)

[H2, w2] = freqz(c,1,512);         % Frequency response of filter
%freqz(c,1,512)

[H3, w3] = freqz(d,1,512);         % Frequency response of filter
%freqz(d,1,512)

[H4, w4] = freqz(e,1,512);         % Frequency response of filter
freqz(e,1,512)

%plot(w1/pi, abs(H1), w2/pi, abs(H2),  f, m, w3/pi, abs(H3), w4/pi, abs(H4), f2, m2)
%legend('N = 34, cutoff=0.3', 'N = 200, cutoff=0.3', 'ideal 1', 'N= 34, cutoff=0.65', 'N=34, cutoff = 0.65, window= Chebyshev', 'ideal 2',  'Location', "EastOutside")
 
output = filtfilt(b,1,y2);       % Zero-phase digital filtering
%figure;                       
%subplot(3,2,1); plot(y2,'b'); title('Original Signal')
%subplot(3,2,2); plot(y,'r'); title('Original Signal with No Noise')
%subplot(3,2,3); plot(output,'g'); title('Filtered Signal (N=34), cutoff = 0.3') 

output2 = filtfilt(c,1,y2);       % Zero-phase digital filtering                      
%subplot(3,2,4); plot(output2,'g'); title('Filtered Signal (N=200, cutoff = 0.3)') 

output3 = filtfilt(d,1,y2);       % Zero-phase digital filtering                       
%subplot(3,2,5); plot(output3,'g'); title('Filtered Signal (N= 34, cutoff = 0.65)') 

output4 = filtfilt(e,1,y2);       % Zero-phase digital filtering                       
%subplot(3,2,6); plot(output4,'g'); title('Filtered Signal (N= 34, cutoff = 0.65, window = Chebyshev )') 






%% FFT

% Ideal HP filter

close all
clear all
clc 

% Read image and calculate FFT
I = imread('driver.jpg');
F = fftshift(fft2(I));
Fa = abs(F);
Flog = log(1+Fa); %because of scale of absolute values have a large range we need to convert to a log scale
Fmax = max(Flog(:)); %need to find max of flog to eventually normailze log values
figure,
subplot(1,2,1); imshow(I); title('Original image');
subplot(1,2,2); imshow(im2uint8(Flog/Fmax),[]); title(' FFT of original image');
% Create ideal HP filter
[rows, columns] = size(I);
[x,y] = meshgrid (-500:499, -315:314); % x columns and y rows
z = sqrt(x.^2 + y.^2);
c = z >15;

% Create butterworth LP filter
bh = hbutterworth(I,15,1);

%Pass image through filter in frequnecy domain

F_smooth = F .* bh;

% Pass image through filter in frequnecy domain

F_filtered = F .* c;

% Restore filtered image
figure,
subplot(2,2,1); imshow(c,[]);

F_inv = ifft2(F_filtered);
F_ia = abs(F_inv);
F_iam = max(F_ia(:));

subplot(2,2,2); imshow(F_ia/F_iam);

F_inv2 = ifft2(F_smooth);
F_ia2 = abs(F_inv2);
F_iam2 = max(F_ia2(:));
subplot(2,2,4); imshow(F_ia2/F_iam2);

subplot(2,2,3);imshow(bh,[]);

function bl = lbutterworth(im,d,n)
	%Create LP butterworth filter which size is same 
	% as the size of the image
	[h,w,f] = size(im);
	[x,y] = meshgrid(-floor(w/2) : floor((w-1)/2), - floor(h/2) : floor((h-1)/2));
	bl = 1./(1+(sqrt(2)-1) * ((x.^2 + y.^2)/d^2).^n);
end

function bh = hbutterworth(im,d,n)
	%Create LP butterworth filter which size is same 
	% as the size of the image
	bh = 1 - lbutterworth(im,d,n);
end
