function spectrum=analyzeFFT(raw,traceName,options)
%ANALYZEFFT Compute a one-sided FFT from an LTspice transient trace.
arguments
 raw; traceName (1,1) string
 options.StartTime_s (1,1) double=0; options.StopTime_s (1,1) double=inf
 options.SampleCount (1,1) double {mustBeInteger,mustBePositive}=4096
 options.Window (1,1) string {mustBeMember(options.Window,["hann","rectangular"])}="hann"
end
if isa(raw,"radia.ltspice.RawRead"),data=raw.Data;elseif isstruct(raw)&&isfield(raw,"waveform"),data=raw.waveform;elseif isstruct(raw),data=raw;else,error("radia:ltspice:FFTInput","raw must be RawRead or RAW/run struct.");end
if data.names(1)~="time",error("radia:ltspice:FFTAnalysis","FFT requires transient time data.");end
j=find(data.names==traceName,1);if isempty(j),error("radia:ltspice:TraceNotFound","Trace not found: %s",traceName);end
t=real(data.values(:,1)); y=data.values(:,j); keep=t>=options.StartTime_s&t<=options.StopTime_s;
t=t(keep);y=y(keep);if numel(t)<2,error("radia:ltspice:FFTSamples","FFT interval has fewer than two samples.");end
uniform=linspace(t(1),t(end),options.SampleCount).'; samples=interp1(t,y,uniform,"linear"); samples=samples-mean(samples);
if options.Window=="hann",window=.5-.5*cos(2*pi*(0:options.SampleCount-1)'/(options.SampleCount-1));else,window=ones(options.SampleCount,1);end
coherentGain=mean(window); transformed=fft(samples.*window); count=floor(options.SampleCount/2)+1;
amplitude=abs(transformed(1:count))/(options.SampleCount*coherentGain);if count>2,amplitude(2:end-1)=2*amplitude(2:end-1);end
dt=uniform(2)-uniform(1); frequency=(0:count-1)'/(options.SampleCount*dt);
spectrum=struct("schema","radia.ltspice.fft.v1","trace",traceName,"frequency_hz",frequency,"amplitude",amplitude,"sample_time_s",dt,"sample_count",options.SampleCount,"window",options.Window);
end
