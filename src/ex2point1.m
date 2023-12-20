% Use 1000 equispaced points on the interval [-1,1].
t = linspace(-2, 2, 1000);

% Sample a smooth function
y = abs(t) .*(2+cos(t)) .* sign(t);
% Try a non-smooth function also:
% y = abs(t) .* exp(t)
%% decontruct signal, question 2.1

% Compute its wavelet transform, four levels deep, using the Daubechies 2
% wavelet.
[c,l] = wavedec(y, 5, 'db2');

% Visualize the coefficients on a logarithmic scale.
% Try to explain what you see! Experiment with other wavelets and, again,
% try to understand the different results.
figure
semilogy(abs(c));
xlabel("$i$",Interpreter="latex");
ylabel("$|c_i|$",Interpreter="latex");

%% recontruct signal, task 2.2

% apply wavelet transform on signal with added noise
rng(42)
epsilon = 1e-1;
noise = epsilon*rand(size(y));
ynoise = y + noise;
[cnoise,lnoise] = wavedec(ynoise, 5, 'db2');

% plot clean and noisy function
figure
plot(t,y)
hold on
plot(t,ynoise)
hold off
xlabel("$t$",Interpreter="latex");
ylabel("function",Interpreter="latex");
legend("$f(t)$","$f(t)+\epsilon \mathcal{U}(0,1)$",Interpreter="latex")

% Find small coefficents and set them to zero.
%T = max(abs(cnoise));
%delta = 10e-4*T;
delta = 1e-1;

% change here to wanted thresholding method
cnoiseInit = cnoise;
[cnoise,I] = Hard_threshold(delta,cnoise);
%[cnoise,I] = Soft_threshold(delta,cnoise);

% How many coefficients did we put to zero:
length(I)
% out of a total of
length(cnoise)
% ratio put to 0
(1-length(I)/length(cnoise))*100
errCoeff = abs(c-cnoise);

% coefficent error between real and thresholded noisy
figure
semilogy(errCoeff);
xlabel("$i$",Interpreter="latex");
ylabel("$|c_i-\hat{c}_i|$",Interpreter="latex");

% Reconstruct the signal
y2 = waverec(cnoise, lnoise, 'db2');
% Plot the error on a logarithmic scale. Experiment with the threshold
% above and see what the effect is
bias = mean(noise); %mean of uniform distrubution
err = abs(y-y2+bias);
errTotal = norm(err)
errTotalnoise = norm(noise)
figure
semilogy(t, err)
hold on
%semilogy(t, noise)
yline(bias,Label="noise",Interpreter="latex")
hold off
xlabel("$t$",Interpreter="latex");
ylabel("Errors",Interpreter="latex");
legend("$|f(t_i)-\hat{f}(t_i)+mean(noise)|$",Interpreter="latex")

%% Question 2.3, finding best threshold for certain parameters
x = linspace(-10,0,101);
deltaList = 10.^x;
errorList = zeros(size(deltaList));
errorMat = zeros(length(t),length(deltaList))';
F = zeros(length(t),length(deltaList))';
errorCoeffList = zeros(size(deltaList));
for i = 1:length(deltaList)
    delta = deltaList(i);
    %[cnoise,I] = Hard_threshold(delta,cnoiseInit);
    [cnoise,I] = Soft_threshold(delta,cnoiseInit);
    
    errorCoeffList(i) = mse(c,cnoise);
    % Reconstruct the signal
    y2 = waverec(cnoise, lnoise, 'db2');
    F(i,:) = y2-bias;
    errorMat(i,:) = abs(y-y2+bias);
    %err = norm(abs(y-y2+bias));
    err = mse(y,y2+bias);
    errorList(i) = err;
end
[BestErr,index] = min(errorList);
[BestErrCoeff,indexCoeff] = min(errorCoeffList);
bestDelta = deltaList(indexCoeff)
BestErrCoeff

% Plot the error on a logarithmic scale. 
figure
semilogy(t, errorMat(indexCoeff,:))
hold on
%semilogy(t, noise)
yline(bias,Label="noise",Interpreter="latex")
hold off
xlabel("$t$",Interpreter="latex");
ylabel("Errors",Interpreter="latex");
legend("$|f(t_i)-\hat{f}(t_i)+mean(noise)|$",Interpreter="latex")

% plot best recontructed function
figure
plot(t,y)
hold on
plot(t,F(indexCoeff,:))
hold off
xlabel("$t$",Interpreter="latex");
ylabel("function",Interpreter="latex");
legend("$f(t)$","$\hat{f}(t)-mean(noise)$",Interpreter="latex")

% function to apply hard thresholding
function [c,I] = Hard_threshold(delta, c)
    I = find(abs(c) < delta);
    c(I) = 0;
end

% function to apply soft thresholding
function [c, I] = Soft_threshold(delta,c)
    I = find(abs(c) < delta);
    c = sign(c).*(abs(c)-delta);
    c(I) = 0;
end