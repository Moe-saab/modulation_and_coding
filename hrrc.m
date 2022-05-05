function [hrrc_time,hrrc_frequency] = hrrc(Fs,Tsymb,Ntaps,rolloff)

	frequency_axis = linspace(-Fs/2,Fs/2,Ntaps); 
	%since every sampling frequency the spectrum repeats, so we plot in frequency domain
	%over 1 sampling frequency enough. Now since we are plotting in frequency domain first
	%we also have to sample in frequency domain, that is what we define as Ntaps
	%so we just take linearly space Ntaps samples in the Fs 
	
	%raised cosine filter in frequency domain
	raised_cosine_filter_f = zeros(1,Ntaps);
    i=1;
    for f = frequency_axis
	    if (abs(f)<(1-rolloff)/(2*Tsymb) )
	       raised_cosine_filter_f(i)=Tsymb;
	    elseif(abs(f)<=(1+rolloff)/(2*Tsymb))
	       raised_cosine_filter_f(i)=Tsymb*(1+cos(pi*Tsymb*(abs(f)-(1-rolloff)/(2*Tsymb))/rolloff))/2;
	    else
	       raised_cosine_filter_f(i)=0;
	    end
	    i=i+1;
    end

    hrrc_frequency = sqrt(raised_cosine_filter_f);
    
    hrrc_time = fftshift(ifft(ifftshift(hrrc_frequency)));
    hrrc_time = hrrc_time/sqrt(sum(hrrc_time.^2)); %normalization of amplitude
    hrrc_time = real(hrrc_time);
    time_axis = 1/Fs * (-Ntaps/2:Ntaps/2 -1);
    
    %plot the filter in time and frequency domain
%     figure;
%     subplot(1,2,1);
%     plot(time_axis,hrrc_time,'k');
%     title('Half root raised cosine filter impulse respones');
%     xlabel('time sec');
%     hold on;grid on;plot(time_axis,hrrc_time,'y*');
%     
%     subplot(1,2,2);
%     plot(frequency_axis,abs(hrrc_frequency),'k');
%     title('Half root raised cosine filter transfer function');
%     xlabel('Hz');
%     hold on;grid on;plot(frequency_axis,hrrc_frequency,'y*');
    
end
