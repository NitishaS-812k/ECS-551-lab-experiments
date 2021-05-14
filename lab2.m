%pcm receiver

%input
Yd = Xd; %without noise

%sipo
Ydpt = reshape(Yd, [n, length(Yd)/n]);  %convert serial bits to parallel
Ydp = Ydpt';

%decoding
IndYt = bi2de(Ydp);
IndY = IndYt';
minL = levels(1);
gapL = levels(2)-levels(1);
Ys = minL + gapL* IndY;      %conversion to levels

%sample and hold
Rat = int64(tPcgrid/tG);      %number of time points between samples
Ysh = zeros(1, length(Ys)*Rat);  %initialise S/H waveform
for k = 1:length(Ys)
    Ysh((k-1)*Rat+1:k*Rat) = Ys(k);
end
Ysh = Ysh(1:(length(Ys)-1)*Rat+1);   %ignore waveform for last sample

%low pass filtering
[b,a] = ellip(2,1,80,0.05);  
Yshf = filter(b,a,Ysh);
%data is sampled at Rat*fs = 100kHz, cut-off frequency = fs/2 = 5kHz
%normalized cutoff frequency = 5/100 = 0.05

%tackle amplitude and delay variation due to imperfect filtering
AdjGain = 1.5;
AdjDelay = -13; 
Yshfp = AdjGain*circshift(Yshf, AdjDelay);

figure(1);
subplot(2,1,1);
stem(tPs(1:5), Xq(1:5), 'b', 'LineWidth', 2);
grid on
xlim([0, 4.5*10^(-4)]);
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Quantized Signal');

subplot(2,1,2);
stem(tPs(1:5), Ys(1:5), 'black', 'LineWidth', 2);
grid on
xlim([0, 4.5*10^(-4)]);
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Decoded Signal');

figure(2);
subplot(2,1,1);
stem(tPs, Ys, 'black', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Decoded Signal');

subplot(2,1,2);
stairs(tPc, Ysh, 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Sample and Hold waveform');

figure(3);
subplot(2,1,1);
plot(tPc, Yshf, 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Reconstructed signal');

subplot(2,1,2);
plot(tPc, Yshfp, 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Shifted signal');