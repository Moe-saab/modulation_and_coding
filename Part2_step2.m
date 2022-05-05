%% Mohamad Saab, Mehmet Fatih, Muhamad Tekin
  %Part2_step2 comparing Hard decoding BER vs no decoding and encoding
%% 
clear all;close all;clc;
%part1
%nbits = 1000;
Nbps = 1;
cutoffreq = 1e6;
rolloff = 0.3;
m = 4; %upsampling ratio
Ntaps =100; %no of taps 
EbN0 = -3:1:15; %snr
Tsymb = 1/(2*cutoffreq);
symbrate = 1/Tsymb;
Fs = symbrate*m;    %sampling frequency
BlockSize = 128;
BlockNb=10*Nbps;
CodeRate = 1/2;
nbits = BlockSize*BlockNb; 
H0 = makeLdpc(BlockSize, BlockSize/CodeRate,0,1,3);  %128*256
repeat = 10;
average_uncoded_BER = zeros(1,length(EbN0));
average_coded_BER = zeros(1,length(EbN0));
for r = 1:repeat

%bit generation
bits_tx = randi(2,1,nbits)-1; %generate sequence of 0's and 1's


%LDPC
blocks=reshape(bits_tx,BlockSize,BlockNb);          
[checkbits, H] = makeParityChk(blocks, H0, 0); % Generate parity check vector bases on LDPC matrix H
                                            %  checkbits       : Check bits 
                                            % H : LDPC matrix
blocks=blocks.';
checkbits=checkbits.';

codedbits=horzcat(checkbits,blocks); %4x256 or BlockNb x 2Blocksize
codedbits_tx=reshape(codedbits.',[],1); %make column vector of size 1024


%mapping
if Nbps>1
    signal_tx_uncoded = mapping(bits_tx.', Nbps, 'qam').';
    signal_tx_coded = mapping(codedbits_tx, Nbps, 'qam').';
else
    signal_tx_uncoded = mapping(bits_tx.', Nbps, 'pam').';
    signal_tx_coded = mapping(codedbits_tx, Nbps, 'pam').';
end

%upsampling
upsampled_uncoded_signal = upsample(signal_tx_uncoded,m);
upsampled_coded_signal = upsample(signal_tx_coded,m);

%hRRC at Tx
% [h_RRC,H_RRC] =  RRC(Fs,Tsymb,N,rolloff,Nbps,1,m);
[hrrc_time,hrrc_frequency] = hrrc(Fs,Tsymb,Ntaps,rolloff);
filtered_uncoded_signal_tx = conv(upsampled_uncoded_signal,hrrc_time);
filtered_coded_signal_tx = conv(upsampled_coded_signal,hrrc_time);


% figure("Name","TX coded vs uncoded_signal");
% subplot(2,3,1);
% plot(signal_tx_uncoded,'g*');
% title("Symblos to be transmitted before upsampling");
% subplot(2,3,2);
% plot(upsampled_uncoded_signal,'b*'); 
% title("upsampled modulated signal","before filtering")
% subplot(2,3,3);
% plot(filtered_uncoded_signal_tx,'r*'); title("filtered signal to be transmited");


% subplot(2,3,4);
% plot(signal_tx_coded,'g*');
% title("Symblos to be transmitted before upsampling");
% subplot(2,3,5);
% plot(upsampled_coded_signal,'b*'); 
% title("upsampled modulated signal","before filtering")
% subplot(2,3,6);
% plot(filtered_coded_signal_tx,'r*'); title("filtered signal to be transmited");




%adding noise
uncoded_signal_energy = trapz((abs(filtered_uncoded_signal_tx).^2))*(1/Fs);
Eb_uncoded = uncoded_signal_energy/(2*nbits);
N0_uncoded = Eb_uncoded./(10.^(EbN0/10));
noise_power_uncoded = 2*N0_uncoded*Fs;
% num of symbols multiplied by the upsampling rate + number of taps
adwg_noise_uncoded = zeros(length(EbN0),length(signal_tx_uncoded)*m+Ntaps-1);
signal_rx_uncoded = zeros(length(EbN0),length(signal_tx_uncoded)*m+Ntaps-1);
adwg_noise_coded = zeros(length(EbN0),length(signal_tx_coded)*m+Ntaps-1);
signal_rx_coded = zeros(length(EbN0),length(signal_tx_coded)*m+Ntaps-1);

for ii = 1:length(EbN0)
adwg_noise_uncoded(ii,:) = sqrt(noise_power_uncoded(ii)/2).*(randn(1,length(signal_tx_uncoded)  * m +Ntaps-1) +1i*randn(1,length(signal_tx_uncoded)*m +Ntaps-1));
signal_rx_uncoded(ii,:) = filtered_uncoded_signal_tx + adwg_noise_uncoded(ii,:);
end

coded_signal_energy = trapz((abs(filtered_coded_signal_tx).^2))*(1/Fs);
Eb_coded = coded_signal_energy/(2*nbits / CodeRate);
N0_coded = Eb_coded./(10.^(EbN0/10));
noise_power_coded = 2*N0_coded*Fs;
% num of symbols multiplied by the upsampling rate + number of taps
for ii = 1:length(EbN0)
adwg_noise_coded(ii,:) = sqrt(noise_power_coded(ii)/2).*(randn(1,length(signal_tx_coded) * m +Ntaps-1) +1i*randn(1,length(signal_tx_coded) * m +Ntaps-1));
signal_rx_coded(ii,:) = filtered_coded_signal_tx + adwg_noise_coded(ii,:);
end



%hRRC at Rx 
a_uncoded_filtered_signal_rx = zeros(length(EbN0),length(signal_tx_uncoded)*m+2*(Ntaps-1));
uncoded_filtered_signal_rx = zeros(length(EbN0),length(signal_tx_uncoded)*m);

a_coded_filtered_signal_rx = zeros(length(EbN0),length(signal_tx_coded)*m+2*(Ntaps-1));
coded_filtered_signal_rx = zeros(length(EbN0),length(signal_tx_coded)*m);

for ii = 1:length(EbN0)
a_uncoded_filtered_signal_rx(ii,:) = conv(signal_rx_uncoded(ii,:), fliplr(hrrc_time));
uncoded_filtered_signal_rx(ii,:) = a_uncoded_filtered_signal_rx(ii,Ntaps:end-(Ntaps-1)); 

a_coded_filtered_signal_rx(ii,:) = conv(signal_rx_coded(ii,:), fliplr(hrrc_time));
coded_filtered_signal_rx(ii,:) = a_coded_filtered_signal_rx(ii,Ntaps:end-(Ntaps-1)); %didn't know why
end


%downsampling
for ii = 1:length(EbN0)
uncoded_downsampled_signal(ii,:) = downsample(uncoded_filtered_signal_rx(ii,:), m);
coded_downsampled_signal(ii,:) = downsample(coded_filtered_signal_rx(ii,:), m);
end

% %figures
% figure("Name","RX coded vs uncoded_signal");
% subplot(2,3,1);
% plot(signal_rx_uncoded,"r*")
% title("uncoded received signal before filtering");
% subplot(2,3,2);
% plot(uncoded_filtered_signal_rx,"y*")
% title("uncoded Filtered RX signal");
% subplot(2,3,3);
% plot(uncoded_downsampled_signal,"g*");
% title("uncoded Downsampled RX signal");
% 
% 
% subplot(2,3,4);
% plot(signal_rx_coded,"r*")
% title("coded received signal before filtering");
% subplot(2,3,5);
% plot(coded_filtered_signal_rx,"y*")
% title("coded Filtered RX signal");
% subplot(2,3,6);
% plot(coded_downsampled_signal,"g*");
% title("coded Downsampled RX signal");



%demapping
uncoded_bits_rx = zeros(length(EbN0),length(bits_tx));
coded_bits_rx = zeros(length(EbN0),length(codedbits_tx));
for ii = 1:length(EbN0)
if Nbps>1
    uncoded_bits_rx(ii,:) = demapping(uncoded_downsampled_signal(ii,:).', Nbps, 'qam');
    coded_bits_rx(ii,:) = demapping(coded_downsampled_signal(ii,:).', Nbps, 'qam');
else
    uncoded_bits_rx(ii,:) = demapping(real(uncoded_downsampled_signal(ii,:).'), Nbps, 'pam');
    coded_bits_rx(ii,:) = demapping(real(coded_downsampled_signal(ii,:).'), Nbps, 'pam');
end
end


%hard_decoding
for ii = 1:length(EbN0)
for j = 1:BlockNb
    codeword = coded_bits_rx(ii,(j-1)*BlockSize/CodeRate+1:j*BlockSize/CodeRate);
    decoded_codeword = hardDecoding(codeword,H,10);
    decoded_bits_rx_HD((j-1)*BlockSize+1:j*BlockSize,ii) = decoded_codeword(BlockSize+1:BlockSize/CodeRate);
end
end


%BER calculation
uncoded_BER = zeros(1,length(EbN0));
coded_BER_HD = zeros(1,length(EbN0));
for ii=1:length(EbN0)
    for j=1:nbits
        if(bits_tx(1,j)~=uncoded_bits_rx(ii,j))
            uncoded_BER(ii) = uncoded_BER(ii) + 1;
        end
%     disp("the error bits without coding are")
%     disp(uncoded_BER);
%     disp("So BER in percent")
%     disp(uncoded_BER/nbits)
        if(bits_tx(1,j)~=decoded_bits_rx_HD(j,ii))
            coded_BER_HD(ii) = coded_BER_HD(ii) + 1;
        end
    end
%     disp("the error bits with coding are")
%     disp(coded_BER_HD);
%     disp("So BER in percent")
%     disp(coded_BER_HD/nbits)


uncoded_BER(ii) = uncoded_BER(ii)/nbits;
coded_BER_HD(ii) = coded_BER_HD(ii)/nbits;

end

average_uncoded_BER = average_uncoded_BER+uncoded_BER;
average_coded_BER = average_coded_BER+coded_BER_HD;
end
average_uncoded_BER = average_uncoded_BER/r;
average_coded_BER = average_coded_BER/r;

figure;
grid on;
semilogy(EbN0,average_uncoded_BER,'r-');
hold on;
semilogy(EbN0,average_coded_BER,'k-');
title('Hard decoding vs no encoding and decoding');
txt = ['Nbps: ' num2str(Nbps)];
text(-2,1e-3,txt)
legend('Uncoded','Hard Decoding with 10 iterations');
xlabel('Eb/N0 [db]'); ylabel('BER');