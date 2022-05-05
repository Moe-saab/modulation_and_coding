%part1
Nbps = 2;
nbits = Nbps * 10000;
cutoffreq = 1e6;
rolloff = 0.3;
fsymb = cutoffreq/Nbps;
Tsymb=1/fsymb;
m = 4; %upsampling ratio
Ntaps =100; %no of taps 
EbN0 = -4:16;
%Tsymb = 1/(2*cutoffreq);
%symbrate = 1/Tsymb;
Fs = fsymb*m;	%sampling frequency


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
figure("Name","TX signal");
subplot(1,3,1);
plot(signal_tx,'g*');
title("Symblos to be transmitted before upsampling");
subplot(1,3,2);
plot(upsampled_signal,'b*'); 
title("upsampled modulated signal","before filtering")
subplot(1,3,3);
plot(filtered_signal_tx,'r*'); title("filtered signal to be transmited");

%adding noise
signal_energy = (trapz(abs(filtered_signal_tx).^2))*(1/Fs);
Eb = signal_energy/(2*nbits);
N0 = Eb./(10.^(EbN0/10));
noise_power = 2*N0*Fs;
% num of symbols multiplied by the upsampling rate + number of taps
%adwg_noise = sqrt(noise_power/2).*(randn(1,nbits/Nbps * m +Ntaps-1) +1i*randn(1,nbits/Nbps * m +Ntaps-1));
%signal_rx = filtered_signal_tx + adwg_noise;

adwg_noise = zeros(length(EbN0),length(signal_tx)*m+Ntaps-1);
signal_rx = zeros(length(EbN0),length(signal_tx)*m+Ntaps-1);

%adwg_noise = sqrt(noise_power/2).*(randn(1,length(filtered_signal_tx)) +1i*randn(1,length(filtered_signal_tx)));
%signal_rx = filtered_signal_tx + adwg_noise;

for j = 1:length(EbN0)
        adwg_noise(j,:) = sqrt(noise_power(j)/2).*(randn(1,length(filtered_signal_tx))+1i*randn(1,length(filtered_signal_tx)));
        signal_rx(j,:) = filtered_signal_tx + adwg_noise(j,:);
end

% adding CFO (frequency offset)
cst = 0;
t1 = 2*pi*cst;
ppm = [0 2 5 10];
cfo = ppm*1e-6*2e9;
ratio = cfo/(Fs);
signal_rx_sync_errors=zeros(length(EbN0),size(signal_rx,2),length(cfo));
for k = 1:length(cfo)
    for i = 1:length(EbN0)
        signal_rx_sync_errors(i,:,k) = signal_rx(i,:).*exp(1j*(2*pi*ratio(k)*(0:length(signal_rx)-1)));
    end
end





% %hRRC at Rx 

filtered_signal_rx_modified = zeros(length(EbN0),length(signal_tx)*m,length(cfo));
filtered_signal_rx_sync_errors = zeros(length(EbN0),length(signal_tx)*m,length(cfo));
for k = 1:length(cfo)
    for i =1:length(EbN0)
        filtered_signal_rx_sync_errors = conv(signal_rx_sync_errors(i,:,k),fliplr(hrrc_time));
        filtered_signal_rx_modified(i,:,k) = filtered_signal_rx_sync_errors(Ntaps:end-(Ntaps-1));
    end                                                                                     %
end 

%downsampling
downsampled_signal_sync_errors = zeros(length(EbN0),length(filtered_signal_rx_modified)/m,length(cfo));
for k = 1:length(cfo)
    for i =1:length(EbN0)
        downsampled_signal_sync_errors(i,:,k) = downsample(filtered_signal_rx_modified(i,:,k), m);
    end
end


%figures
% figure("Name","RX signal");
% subplot(1,3,1);
% plot(signal_rx,"r*")
% title("received signal before filtering");
% subplot(1,3,2);
% plot(filtered_signal_rx,"y*")
% title("Filtered RX signal");
% subplot(1,3,3);
% plot(downsampled_signal,"g*");
% title("Downsampled RX signal");
% hold on; plot(downsampled_signal_sync_errors,"r*");
% legend('ideal with 8db snr','constellation after adding cfo 2ppm');
% 
% figure;plot(downsampled_signal,"g*")
% title("Downsampled RX signal");
% hold on; plot(downsampled_signal_sync_errors,"r*");
% legend('ideal with 8db snr','constellation after adding cfo 2ppm');
%demapping
demodulated_bits = zeros(length(EbN0),nbits,length(cfo));
    if Nbps>1
        for k = 1:length(cfo)
            for i =1:length(EbN0)
                demodulated_bits(i,:,k) = demapping(downsampled_signal_sync_errors(i,:,k).',Nbps,'qam').';
            end                                                                                     %
        end
    else
        for k = 1:length(cfo)
            for i =1:length(EbN0)
                demodulated_bits(i,:,k) = demapping(downsampled_signal(i,:,k).',Nbps,'pam').';
            end                                                                                     %
        end
    end


%BER calculation
% BER = 0;
% for i=1:nbits
% 	if(bits_tx(i)~=bits_rx(i))
% 		BER = BER + 1;
% 	end
% end
% disp("the error bits are")
% disp(BER);
% disp("So BER in percent")
% disp(BER/nbits)

error_symbols = zeros(length(EbN0),length(cfo));
for k = 1:length(cfo)
    for i =1:length(EbN0)
        error_symbols(i,k) = sum(bits_tx ~= demodulated_bits(i,:,k));
    end                                                                                     %
end

errors(:,:,1) = error_symbols;

ber_cfo = errors/nbits;

figure;
semilogy(EbN0,ber_cfo(:,1));
hold on;
semilogy(EbN0,ber_cfo(:,2),'-');
semilogy(EbN0,ber_cfo(:,3),'-');
semilogy(EbN0,ber_cfo(:,4),'-');
legend('1','2','3','4')













