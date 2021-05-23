%pcm receiver

%input
Yd = Xd; %without noise
Yd_n10 = decision_device(Xd, 10);   %Xd with channel snr 10dB
Yd_n01 = decision_device(Xd, -10);  %Xd with channel snr -10dB

%without channel noise
Ydp = sipo(Yd,n);                       %sipo
Ys = decoding(Ydp,levels);              %decoding
Ysh = sample_and_hold(tPcgrid,tG,Ys);   %sample and hold
Yshf_ell = filtering(Ysh, 'elliptic');  %elliptic filtering
Yshf_butter = filtering(Ysh, 'butter'); %butterworth filtering
Yshfp_ell = adjust_gain_delay(Yshf_ell, 'elliptic'); %elliptic delay and gain adjustment
Yshfp_butter = adjust_gain_delay(Yshf_butter, 'butter'); %butterworth filter delay and gain adjustment

%with 10dB awgn channel noise
Ydp_n10 = sipo(Yd_n10,n);                      %sipo
Ys_n10 = decoding(Ydp_n10,levels);             %decoding
Ysh_n10 = sample_and_hold(tPcgrid, tG,Ys_n10); %sample and hold
Yshfn10_ell = filtering(Ysh_n10, 'elliptic');  %elliptic filtering
Yshfn10_butter = filtering(Ysh_n10, 'butter'); %butter filtering
Yshfpn10_ell = adjust_gain_delay(Yshfn10_ell, 'elliptic'); %elliptic delay and gain adjustment
Yshfpn10_butter = adjust_gain_delay(Yshfn10_butter, 'butter'); %butterworth filter delay and gain adjustment

%with -10dB awgn channel noise
Ydp_n01 = sipo(Yd_n01,n);                       %sipo
Ys_n01 = decoding(Ydp_n01,levels);              %decoding
Ysh_n01 = sample_and_hold(tPcgrid, tG,Ys_n01);  %sample and hold
Yshfn01_ell = filtering(Ysh_n01, 'elliptic');   %elliptic filtering
Yshfn01_butter = filtering(Ysh_n01, 'butter');  %butterworth filtering
Yshfpn01_ell = adjust_gain_delay(Yshfn01_ell, 'elliptic'); %elliptic delay and gain adjustment
Yshfpn01_butter = adjust_gain_delay(Yshfn01_butter, 'butter'); %butterworth filter delay and gain adjustment

%data is sampled at Rat*fs = 100kHz, cut-off frequency = fs/2 = 5kHz
%normalized cutoff frequency = 5/100 = 0.05

%plots for without noise case
figure(1);
%decoded signal plot
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
title('without noise')


figure(2);
%sample and hold waveform plot no channel-noise
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
%sample and hold waveforms for 10dB and -10dB snr 
subplot(2,1,1);
stairs(tPc, Ysh_n10, 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Sample and Hold waveform');
title('with 10dB snr')

subplot(2,1,2);
stairs(tPc, real(Ysh_n01), 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Sample and Hold waveform');
title('with -10dB snr')


figure(4);
%reconstruction plot for no channel noise elliptic filter
subplot(2,1,1);
plot(tPc, Yshf_ell, 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Reconstructed signal');
title('signal reconstructed with elliptical filter(no noise)')

subplot(2,1,2);
plot(tPc, Yshfp_ell, 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Shifted signal');
title('adjusted for elliptical filter(no noise)')



figure(5);
%reconstruction plot for no-noise butterworth filter
subplot(2,1,1);
plot(tPc, Yshf_butter, 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Reconstructed signal');
title('signal reconstructed with butterworth filter(no noise)')

subplot(2,1,2);
plot(tPc, Yshfp_butter, 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Shifted signal');
title('adjusted for butterworth filter(no noise)')


figure(6)
%reconstruction plot for 10dB snr elliptic filter
subplot(2,1,1);
plot(tPc, real(Yshfn10_ell), 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Reconstructed signal');
title('reconstructed signal with elliptic filter(snr = 10dB)');

subplot(2,1,2);
plot(tPc, real(Yshfpn10_ell), 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Shifted signal');
title('adjusted waveform for elliptic filter(snr = 10dB)');


figure(7);
%reconstruction plot for 10dB snr butterworth filter
subplot(2,1,1);
plot(tPc, real(Yshfn10_butter), 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Reconstructed signal');
title('reconstructed signal with butterworth filter(snr = 10dB)');

subplot(2,1,2);
plot(tPc, real(Yshfpn10_butter), 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Shifted signal');
title('adjusted waveform for butterworth filter(snr = 10dB)');

figure(8);
%reconstruction plot for -10dB snr elliptic filter
subplot(2,1,1);
plot(tPc, real(Yshfn01_ell), 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Reconstructed signal');
title('reconstructed signal with elliptic filter(snr = -10dB)');

subplot(2,1,2);
plot(tPc, real(Yshfpn01_ell), 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Shifted signal');
title('adjusted waveform for elliptic filter(snr = -10dB)');


figure(9);
%reconstruction plot for -10dB snr butterworth filter
subplot(2,1,1);
plot(tPc, real(Yshfn01_butter), 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Reconstructed signal');
title('reconstructed signal with butterworth filter(snr = -10dB)');

subplot(2,1,2);
plot(tPc, real(Yshfpn01_butter), 'b', 'LineWidth', 2);
grid on
hold on
plot(tPc, X, 'r', 'LineWidth', 1, 'LineStyle', '--');
hold off
xt = get(gca, 'XTick');
set(gca, 'XTick',xt, 'XTickLabel', xt*1000);
xlabel('Time(ms)');
ylabel('Shifted signal');
title('adjusted waveform for butterworth filter(snr = -10dB)');


function Ys = decoding(Ydp,levels)
    %decoding
    IndYt = bi2de(Ydp);
    IndY = IndYt';
    minL = levels(1);
    gapL = levels(2)-levels(1);
    Ys = minL + gapL* IndY;  %conversion to levels
end

function Ys = decision_device(Yd, snrdB)
    %takes snr in dB to simulate a gaussian random variable to act as
    %channel noise
    Ys = zeros(size(Yd));
    snr = 10^(snrdB/10);
    No = 1/snr;  %assuming Eb = 1
    for i = 1:length(Yd)
        N = sqrt(No/2)*randn(1);
        if Yd(i) + N < 0.5
            Ys(i) = 0;
        end
        if Yd(i) + N >= 0.5
            Ys(i) = 1;
        end
    end   
end

function Ysh = sample_and_hold(tPcgrid,tG,Ys)
    %sample and hold
    Rat = int64(tPcgrid/tG);         %number of time points between samples
    Ysh = zeros(1, length(Ys)*Rat);  %initialise S/H waveform
    for k = 1:length(Ys)
        Ysh((k-1)*Rat+1:k*Rat) = Ys(k);
    end
    Ysh = Ysh(1:(length(Ys)-1)*Rat+1);   %ignore waveform for last sample
end

function Ydp = sipo(Yd,n)
    Ydpt = reshape(Yd, [n, length(Yd)/n]);  %convert serial bits to parallel
    Ydp = Ydpt';
end

function Yshf = filtering(Ysh, filter_type)
    %low pass filtering
    if strcmp(filter_type,'elliptic')  %elliptic filtering
        [b,a] = ellip(2,1,80,0.05);  
        Yshf = filter(b,a,Ysh);
    end
    if strcmp(filter_type, 'butter')   %butterworth filtering
       [b,a] = butter(3,0.05);  
       Yshf = filter(b,a,Ysh);
    end
end 

function Yshfp = adjust_gain_delay(Yshf, filter_type)
    %tackle amplitude and delay variation due to imperfect filtering
    if strcmp(filter_type,'elliptic')
        AdjGain = 1.5;
        AdjDelay = -13; 
        Yshfp = AdjGain*circshift(Yshf, AdjDelay);
    end
    if strcmp(filter_type, 'butter')
        AdjGain = 2;
        AdjDelay = 20;
        Yshfp = AdjGain*circshift(Yshf, AdjDelay);
    end
end