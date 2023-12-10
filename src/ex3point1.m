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
Ared(mask) = 0;
Ared_guess = Ared;
Ared_extr = [Ared(49,70:80) ; Ared(61,70:80)];
Ared_extr(:,2:12) = interp1([1 4],Ared_extr(:,[1 13]).',[2 3]).';
n= size(Ared_extr,1);
data_interp = interp1(1:n,Ared_extr,linspace(1,n,12*(n-1)+1),'spline');
Ared_guess(50:60,70:80)=data_interp(2:12,:);

%scheme
%Ared + ()
