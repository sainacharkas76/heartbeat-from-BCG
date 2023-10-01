function sigfilt=filtHP(signal,fs)

% high pass - TO remove respiration "noise"

    [b,a]=butter(4,5/(fs/2),'high');
    sigfilt=filtfilt(b,a,signal);
    [b,a]=butter(4,25/(fs/2),'low');
    sigfilt=filtfilt(b,a,sigfilt);
    
end
