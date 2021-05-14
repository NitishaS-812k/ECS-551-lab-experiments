close all;
clc;
num_bit = 10000;    %signal length
max_run = 20;       %max number of iterations for a single snr
Eb = 1;             %bit energy
SNRdB = -10:2:10;   %signal to noise ratio(in dB)
SNR = 10.^(SNRdB/10); %convert snr to normal scale
BER_sim = zeros(size(SNR));
BER_the = zeros(size(SNR));

for count = 1:length(SNR)    %beginning of loop for different snr
    avgError = 0;
    for run_time = 1:max_run  %beginning of loop for different runs
        Error = 0;
        D = randi([0 1], 1, num_bit); %baseband modulation
        X = 2*D - 1;
        No = Eb/SNR(count);
        N = sqrt(No/2)*randn(1, num_bit);
        Y = X + N;
        for n = 1:num_bit
            if ((Y(n) > 0 && D(n) == 0) ||(Y(n) < 0 && D(n) == 1))
                Error = Error + 1;
            end
        end
        Error = Error/num_bit;     %calculate error/bit for different runs
        avgError = avgError + Error;    
    end
    BER_sim(count) = avgError/max_run;    %simulated BER
    BER_the(count) = 0.5*erfc(sqrt(SNR(count))); %theoretical BER
end

figure(1);
semilogy(SNRdB, BER_the, 'r', 'LineWidth', 2);
grid on
hold on
semilogy(SNRdB, BER_sim,'marker', '.','color','k', 'MarkerSize', 20);
hold off
legend({'theoretical', 'simulated'});
ylabel('BER');
xlabel('SNR(in dB)');        