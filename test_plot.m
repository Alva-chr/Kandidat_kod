% simy_data_plot.m
% ----------------
% Define your data as [frequency, dB] pairs
data = [
    1   -188.764658;
    2   -191.717359;
    4   -174.187139;
    6   -172.914656;
    8   -183.891174;
   10   -177.794865;
   12   -160.199172;
   14   -159.386479;
   16   -157.747956;
   18   -151.915263;
   20   -117.421005;
   22   -175.219153;
   24   -139.479597;
   26   -124.002068;
   28   -121.698759;
   30   -130.507753;
    3   -185.588023;
    5   -169.295298;
    7   -176.831589;
    9   -180.134006;
   11   -168.944957;
   13   -157.547835;
   15   -165.441883;
   17   -134.311064;
   19   -163.792024;
   21   -160.135861;
   23   -113.687681;
   25   -130.707274;
   27   -132.819582;
   29   -133.218141
];

% Separate into vectors
freq = data(:,1);
db   = data(:,2);

% Sort by frequency (just in case)
[freq, idx] = sort(freq);
db = db(idx);

% Create the plot
figure;
plot(freq, db, '-o', 'LineWidth', 1.5, 'MarkerSize', 6);
xlabel('Frequency');
ylabel('Amplitude (dB)');
title('simy Data vs. Frequency');
grid on;
