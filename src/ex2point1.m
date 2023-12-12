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
ynoise = y + epsilon*rand(size(y));
[cnoise,lnoise] = wavedec(ynoise, 5, 'db2');

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

% change here what thresholding method wanted
[cnoise,I] = Hard_threshold(delta,cnoise);
%[cnoise,I] = Soft_threshold(delta,cnoise);

% How many coefficients did we put to zero:
length(I)
% out of a total of
length(cnoise)
% ratio
length(I)/length(cnoise)
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
err = abs(y-y2);
noise = mean(abs(y-ynoise));
errTotal = norm(err)
errTotalnoise = norm(abs(y-ynoise))
figure
semilogy(t, err)
hold on
%semilogy(t, noise)
yline(noise,Label="noise",Interpreter="latex")
hold off
xlabel("$t$",Interpreter="latex");
ylabel("Errors",Interpreter="latex");
legend("$|f(t_i)-\hat{f}(t_i)|$",Interpreter="latex")

%% Question 2.3

deltaList = linspace(1e-3,1e-10);

for i = 1:length(deltaList)
    
end

function [c,I] = Hard_threshold(delta, c)
    I = find(abs(c) < delta);
    c(I) = 0;
end

function [c, I] = Soft_threshold(delta,c)
    I = find(abs(c) < delta);
    c = sign(c).*(abs(c)-delta);
    c(I) = 0;
end