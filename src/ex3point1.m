% The colormap is empty.
[A,cmap] = imread('bib.png');
A=im2double(A);
figure
imshow(A)

% The tensor contains red, green and blue components
Ared = A(:,:,1);
Agreen = A(:,:,2);
Ablue = A(:,:,3);
Ared_orig= Ared;


mask = zeros(size(Ared));
mask(50:60,70:80) = 1;
mask = mask > 0;
Ared(mask) = 0;  % masked
% Ared(mask_compl) = 0;
Ared_guess = Ared;%B0
Ared_extr = [Ared(49,70:80) ; Ared(61,70:80)];
n= size(Ared_extr,1);
Ared_interp = interp1(1:n,Ared_extr,linspace(1,n,12*(n-1)+1),'spline');
Ared_guess(50:60,70:80)=Ared_interp(2:12,:); %B0

%% one time step of the interation scheme 
% take the first guess Ared_guess, deconstruct and
% reconstruct, take the maskcompliment and add the 
% masked Ared 

[c1,l1] = wavedec2(Ared_guess, 4, 'bior4.4');
[c2,l2] = wavedec2(Agreen, 4, 'bior4.4');
[c3,l3] = wavedec2(Ablue, 4, 'bior4.4');

% Compress by putting small coefficients to zero.
% Experiment with the threshold p.
T = max([max(abs(c1)), max(abs(c2)), max(abs(c3))]);
p = 1e-2;
threshold = p*T
I1 = find(abs(c1) < threshold);
c1(I1) = 0;
I2 = find(abs(c2) < threshold);
c2(I2) = 0;
I3 = find(abs(c3) < threshold);
c3(I3) = 0;

A2red = waverec2(c1, l1, 'bior4.4');
A2green = waverec2(c2, l1, 'bior4.4');
A2blue = waverec2(c3, l1, 'bior4.4');


% A2red should masked complimentarily (I-P)

A2red(mask_compl)=0;

B1 = Ared + A2red; % Ared is already masked 
