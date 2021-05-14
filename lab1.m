%pcm transmitter 

%input
tG = 10^(-5);     %gap between two time points
mint = 0;         %minimum of time scale
maxt = 2*10^(-3); %maximum of time scale
tPc = mint:tG:maxt; %time scale for plot
X = cos(2*pi*10^3*tPc) + 3*sin(2*pi*3*10^3*tPc);

%sampling
fs = 10*10^3;     %sampling frequency
tPcgrid = 1/fs;    %sampling interval
tPs = mint:tPcgrid:maxt; %sampled time axis
tpen = mint:tPcgrid:maxt*3;
Xs = cos(2*pi*10^3*tPs) + 3*sin(2*pi*3*10^3*tPs);   %sampled signal

%quantization
bounds = -3:1:3;    %range boundaries for quantization
levels = -3.5:1:3.5;  %levels of quantization
[ind, Xq] = quantiz(Xs,bounds,levels); 

%encoding
M = length(levels);  %number of quantization levels
n = log2(M);         %No. of bits per sample
Xdp = de2bi(ind,n);  %encoder output, parallel bits

%PISO
lXd = numel(Xdp);    %length of the binary output
Xd = reshape(Xdp', [1, lXd]);  %convert parallel bits to serial

figure(1);
 
subplot(2,1,1);
stem(tPs, Xs, 'b', 'LineWidth', 2);
hold on 
grid on
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xlabel('Time(ms)');
ylabel('Sampled Signal');

subplot(2,1,2);
stem(tPs, Xq,'black', 'LineWidth', 2);
hold on
grid on
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xlabel('Time(ms)');
ylabel('Quantized Signal');

figure(2);
subplot(2,1,1);
stem(tPs(1:5), Xq(1:5), 'b', 'LineWidth', 2); 
grid on
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlim([0, 4.5*10^(-4)]);
xlabel('Time(ms)');
ylabel('Quantized signal');

subplot(2,1,2);
stairs(tpen(1:15),Xd(1:15), 'black', 'LineWidth', 2);
grid on
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
ylim([-0.5, 1.5]);
xlabel('Time(ms)');
ylabel('Bits Transmitted');