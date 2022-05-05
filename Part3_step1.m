%% Author Mohamad Saab
%% Impact of CFO on bit error rate (only ISI)
% Initialization
nbits = 1000;
Nbps = 2;
cutoffreq = 1e6;
rolloff = 0.3;
m = 4; %upsampling ratio
Ntaps =100; %no of taps 
EbN0 = -4:16; %snr
Tsymb = 1/(2*cutoffreq);
symbrate = 1/Tsymb;
Fs = symbrate*m;	%sampling frequency
repeatnumber = 10;
for repeat = 1:repeatnumber
    
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
        adwg_noise(j,:) = sqrt(NoisePower(j)/2).*(randn(1,length(signal_tx)*m+Ntaps-1)+1i*randn(1,length(signal_tx)*m+Ntaps-1));
        signal_rx(j,:) = filtered_signal_tx + adwg_noise(j,:);
    end

end
