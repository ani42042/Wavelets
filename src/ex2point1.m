% Use 1000 equispaced points on the interval [-1,1].
t = linspace(-2, 2, 1000);

% Sample a smooth function
y = abs(t) .*(2+cos(t)) .* sign(t);
% Try a non-smooth function also:
% y = abs(t) .* exp(t)

%% decontruct signal

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

%% recontruct signal

% Find small coefficents and set them to zero.
T = max(abs(c));
I = find(abs(c) < 1e-4*T);
c(I) = 0;

% How many coefficients did we put to zero:
length(I)
% out of a total of
length(c)

% Reconstruct the signal
y2 = waverec(c, l, 'db2');

% Plot the error on a logarithmic scale. Experiment with the threshold
% above and see what the effect is
figure
semilogy(t, abs(y-y2))