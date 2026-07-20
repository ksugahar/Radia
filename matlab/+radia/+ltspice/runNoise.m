function result=runNoise(netlistFile,options)
%RUNNOISE Run an LTspice .noise analysis and expose noise spectra.
arguments
 netlistFile (1,1) string {mustBeFile}; options.Parameters (1,1) struct=struct(); options.Executable (1,1) string=""; options.OutputDirectory (1,1) string=""
end
result=radia.ltspice.run(netlistFile,Parameters=options.Parameters,Executable=options.Executable,OutputDirectory=options.OutputDirectory,RawFormat="binary");
names=result.waveform.names; frequency=real(result.waveform.values(:,1)); noise=struct();
for k=2:numel(names)
 key=matlab.lang.makeValidName(names(k),"ReplacementStyle","hex"); noise.(key)=result.waveform.values(:,k);
end
result.schema="radia.ltspice.noise.v1"; result.frequency_hz=frequency; result.noise_traces=noise;
end
