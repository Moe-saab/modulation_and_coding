%% Author Mohamad Saab
%% Impact of CFO on bit error rate (only ISI)
clear all;
close all;
clc;
% Initialization
nbits = 100000;
Nbps = 4;
cutoffreq = 1e6;
fc = 2e9; %carrier frequency
rolloff = 0.3;
m = 4; %upsampling ratio
Fs = 2*cutoffreq*m; %sampling frequency
Ts = 1/Fs; % sampling time
fsymb = 2*cutoffreq; % symbol frequency
Tsymb = 1/fsymb; % time between two symbols
Ntaps = 31;  %no of taps 
EbN0 = -4:16; %snr

repeatnumber = 3;

ppm = linspace(0,10,6); % we will simulate up to 10 ppm

CFO = ppm*fc*1e-6;

errors = zeros (length(EbN0),length(CFO),repeatnumber);

for t = 1:repeatnumber
    
    %bit generation
    bits_tx = randi(2,1,nbits)-1; %generate sequence of 0's and 1's
    
    %mapping
    if Nbps>1
	signal_tx = mapping(bits_tx.', Nbps, 'qam').';
    else
	signal_tx = mapping(bits_tx.', Nbps, 'pam').';
    end
    
    %upsampling
    upsampled_signal = upsample(signal_tx,m);
    
    %hRRC at Tx
    [hrrc_time,hrrc_frequency] = hrrc(Fs,Tsymb,Ntaps,rolloff);
    filtered_signal_tx = conv(upsampled_signal,hrrc_time);
    
    %adding noise
    signal_energy = trapz(abs(filtered_signal_tx).^2)*(1/Fs);
    Eb = signal_energy/(2*nbits);
    N0 = Eb./(10.^(EbN0/10));
    noise_power = 2*N0*Fs;
    
    adwg_noise = zeros(length(EbN0),length(signal_tx)*m+Ntaps-1);
    signal_rx = zeros(length(EbN0),length(signal_tx)*m+Ntaps-1);

    for j = 1:length(EbN0)
        adwg_noise(j,:) = sqrt(noise_power(j)/2).*(randn(1,length(signal_tx)*m+Ntaps-1)+1i*randn(1,length(signal_tx)*m+Ntaps-1));
        signal_rx(j,:) = filtered_signal_tx + adwg_noise(j,:);
    end
    
    
    % adding CFO (frequency offset)
    t1 = ((0:length(signal_rx)-1))*Ts;
    signal_rx_sync_errors=zeros(length(EbN0),size(signal_rx,2),length(CFO));
    
    for k=1:length(CFO)
        for i = 1:length(EbN0)
            signal_rx_sync_errors(i,:,k) = signal_rx(i,:).*exp(1j*(2*pi*CFO(k)*t1));
        end
    end
 
    
    %hRRC at Rx 
    filtered_signal_rx = zeros(length(EbN0),length(signal_tx)*m,length(CFO));
    filtered_signal_rx_modified = zeros(length(EbN0),length(signal_tx)*m,length(CFO));

    for k = 1:length(CFO)
        for i =1:length(EbN0)
            filtered_signal_rx = conv(signal_rx_sync_errors(i,:,k),fliplr(hrrc_time));
            filtered_signal_rx_modified(i,:,k) = filtered_signal_rx(Ntaps:end-(Ntaps-1));
        end                                                                                     %
    end 

    % Compensation of CFO (in order to observe ISI effect only)
    t2 = ((0:length(filtered_signal_rx_modified)-1))*Ts;
    for k=1:length(CFO)
        for i = 1:length(EbN0)
            filtered_signal_rx_modified(i,:,k) = filtered_signal_rx_modified(i,:,k).*exp(-1j*(2*pi*CFO(k)*t2));
        end
    end

    % Downsampling
    downsampled_signal = zeros(length(EbN0),length(filtered_signal_rx_modified)/m,length(CFO));
    for k = 1:length(CFO)
        for i =1:length(EbN0)
            filtered_signal_rx = conv(signal_rx_sync_errors(i,:,k),fliplr(hrrc_time));
            downsampled_signal(i,:,k) = downsample(filtered_signal_rx_modified(i,:,k),m);
        end                                                                                     %
    end

   

    % Demapping
    demodulated_bits = zeros(length(EbN0),nbits,length(CFO));
    if Nbps>1
	    for k = 1:length(CFO)
            for i =1:length(EbN0)
                demodulated_bits(i,:,k) = demapping(downsampled_signal(i,:,k).',Nbps,'qam').';
            end                                                                                     %
        end
    else
        for k = 1:length(CFO)
            for i =1:length(EbN0)
                demodulated_bits(i,:,k) = demapping(downsampled_signal(i,:,k).',Nbps,'pam').';
            end                                                                                     %
        end
    end
    
    % Finding The Probability of Error in Bits
    error_symbols = zeros(length(EbN0),length(CFO));

    for k = 1:length(CFO)
        for i =1:length(EbN0)
            error_symbols(i,k) = sum(bits_tx ~= demodulated_bits(i,:,k));
        end                                                                                     %
    end



    errors(:,:,t) = error_symbols;

    
end

 ber_qam_4 = (mean(errors,3))/nbits;

 ber_qam_4_theoretical = berawgn(EbN0,'qam',16); % 4-QAM BER Curve


figure;
semilogy(EbN0,ber_qam_4(:,1),'-');
hold on
semilogy(EbN0,ber_qam_4(:,2),'-');
semilogy(EbN0,ber_qam_4(:,3),'-');
semilogy(EbN0,ber_qam_4(:,4),'-');
semilogy(EbN0,ber_qam_4(:,5),'-');
semilogy(EbN0,ber_qam_4(:,6),'-');
semilogy(EbN0,ber_qam_4_theoretical,'o')
legend('1','2','3','4','5','6','7')
ylim([1e-5 1]);