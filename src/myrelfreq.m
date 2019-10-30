function [freq,cfreq] = myrelfreq (x,xbins)

% Return the frequency (FREQ) and cumulative frequency (CFREQ) of the
% values in X over the intervals defined in XBINS
%
% [FREQ,CFREQ] = myrelfreq(X,XBINS)
%
% Sebastian Jaramillo-Riveri
% June, 2013
    
freq = histc(x,xbins);
freq = (freq)./(sum(freq));
cfreq = cumsum(freq);

end