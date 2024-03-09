Now, remove some lines of the waveform*.xml: 
<waveformArm>
....
</waveformArm>

we need to add the automated processing Framework to the waveform*.xml:
    <waveformVectorArm>
        <phaseRequest>
            <beginPhase>ttp</beginPhase>
            <beginOffset>
                <unit>MINUTE</unit>
                <value>-2</value>
            </beginOffset>
            <endPhase>ttp</endPhase>
            <endOffset>
                <unit>MINUTE</unit>
                <value>3</value>
            </endOffset>
        </phaseRequest>
        <fdsnDataSelect>
            <host>service.ncedc.org</host>
        </fdsnDataSelect>
	<bestChannelAtStation/>
        <fullCoverage/>
        <printlineSeismogramProcess/>
        <sacWriter/>
        <responseGain/>
        <rMean/>
        <rTrend/>
        <integrate/>
                <sampleSyncronize/>
                <vectorTrim/>
                <rotateGCP/>
        <sacWriter>
            <workingDir>processedSeismograms</workingDir>
        </sacWriter>
        <legacyExecute>
            <command>echo Sod saved this file</command>
        </legacyExecute>

    </waveformVectorArm>

