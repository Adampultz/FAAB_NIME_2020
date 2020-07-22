
s = Server.default;

s.options.numAnalogInChannels = 8;
s.options.memSize_(8192*16);
s.options.maxLogins = 4;
s.options.headphoneLevel = -6; // sets the headphone level to -6 dB
s.options.speakerMuted = 0;    // enable the speaker amp
s.options.numMultiplexChannels = 0; // do not enable multiplexer channels
s.options.blockSize = 32;
s.options.numInputBusChannels = 4;
s.options.numOutputBusChannels = 2;

s.options.postln;

~numIns = 4;
~anIns = [4,7,5,6];
~stringAmps = [4,2,1,1];
// ~stringAmps = [1];
~syncAmps = Bus.control(s,4); // Amp for scaling incoming amplitudes
~syncFreq = Bus.control(s,4); // Frequency control for oscillators
~variance = Bus.control(s,1).set(0.0001); // Set initial variance for sync
~weights = Bus.control(s,8).setn(8.collect({rrand(0.5, 1.5)})); // Set initial weights for sync
~phaseFreq = Bus.control(s,1).set(0.0001); // Set initial frequency of phasor
~trig = Bus.control(s,1).set(1); 
~distRange = Bus.control(s,1); // Distrubtion range for sum of amplitudes~distRange = Bus.control(s,1).set(10); // Distrubtion range for sum of amplitudes
~countReset = Bus.control(s,1).set(0);
~phasorGate = Bus.control(s,1).set(1);
d = Dictionary.new; // Make a dictionary for holding all information about the GA
4.collect({|i| d.add((\amp_ ++ i).asSymbol -> [0])});
d.add(\mainAmp -> [0]);
d.add(\varSum -> [0]);
d.add(\ampSync -> Array.fill(4,{0}));
d.add(\ampSyncMono -> Array.fill(8,{0}));
~searchSpace = 8.collect({|i| i}); // Make the search space (8 individual weights)
~tempSS = Array.fill(8, 0); // Temporary search space for storing parent indexes
8.collect({|i| d.add((\parent_ ++ i).asSymbol -> rand(100000).asBinaryDigits(17))}); // Generate 8 random weights of 17 binary digits each
d.add(\counter -> [0]);
d.add(\variance -> 0);
d.add(\fitnessLS -> Array.fill(8,{rand(1.0)}).normalizeSum);


o = OSCFunc({|msg|
var mean = msg[3..6].mean; // Amplitude mean of strings
var variance = (4.collect({|i| ((msg[i+3]) - (mean)).pow(2)}).sum/4); // Variance
d[\mainAmp].put(0, d[\mainAmp][0] + mean); // Accumulate mean amplitude of all strings
d[\varSum].put(0, d[\varSum][0] + variance); // Accumulate 

4.collect({|i| d[(\amp_ ++ i).asSymbol].put(0,d[(\amp_ ++ i).asSymbol][0]+ msg[3+i])}); // Running sync to sum of amplitudes
	d[\counter].put(0, msg[9].asInt);
		if(msg[7] < 0.9,
		{~trig.set(0);
		~phasorGate.set(1); // Keep phasor gate open
			~countReset.set(0)},
		{
		~phasorGate.set(-1); // Keep phasor at 0.91
			if(variance > 0.1, // If distribution range is large
				{
				~exec.value; // Execute script for new weights and phasor freq
			})
		}
	)
},'/phasor');

s.waitForBoot {
	"Server Booted".postln;
	(
	SynthDef(\strings,
			{ |rmsLP = 4, filFrq =0.1, lim = 3.0, oscMax =0.15, ampSc =4, ampDiffScale =1, dcScale = 1,  cD=0.2, cDecay= 0.2, dynLpf =8000|
		var locIn = ~syncAmps.kr;
		var sr = SampleRate.ir;
				var anIn = A2K.kr(AnalogIn.ar(~anIns));
				var ctlSr = ControlRate.ir;
				var anScale = LinLin.kr(anIn.round(0.001), 0, 0.827, -60.0, 0.0).dbamp; // Scaled to dB and then linear amplitude
				var pickIn= SoundIn.ar([0,1,2,3]) * anScale * ~stringAmps * 10;
			var sIn = XFade2.ar(pickIn, (pickIn * locIn) * 0.5, A2K.kr(AnalogIn.ar(DC.ar(0))).linlin(0.004, 3.3 / 4.096, -1.0,1.0));
			var sInKr = A2K.kr(pickIn);
				var freq = Pitch.kr(Mix(sInKr))[0]; // Pitch detector
				var freqZero = Pitch.kr(sInKr).clip(30, 4000); // Zero crossings (for crude low-CPU individual pitch estimation)
				var rmsFreq = RMS.kr(freq, 10); // Fast RMS of frequency. Weird but works
				var rmsFreqSlow = RMS.kr(freq, 0.2); // Slow frequency RMS
				var freqDiff = RMS.kr(rmsFreq - rmsFreqSlow, 1).abs; // Absolute difference between fast and slow frequency detection
			var ab = RMS.kr(Mix(pickIn * ampSc), 0.2); // Slow RMS amplitude follower, sum of all inputs
			var rmsQuick = RMS.kr(Mix(pickIn * ampSc), 10); // Fast RMS amplitude follower
				var ampDiff = RMS.kr((ab - rmsQuick), 4).abs; // Absolute difference between slow and fast amplitude RMS
		var a = RMS.kr(pickIn * ampSc, ampDiff.sqrt); // RMS amplitude follower. Amplitude scaled up. Window = Square root of amp difference + 1
		var aL = RMS.kr(SoundIn.ar([0,1,2,3]), ampDiff.reciprocal);
	//	var sync = 4.collect({ |i| Phasor.kr(1, (aL[i].sqrt - [aL[(i+1).mod(4)].sqrt * ~weights.kr(1,(i+1).mod(4)), aL[(i-1).mod(4)].sqrt * ~weights.kr(1,(i-1).mod(4))].mean), 0,ctlSr).fold2(ctlSr /4) / (ctlSr /4) * 1.0 + 1.0});
	var sync = 4.collect({ |i| LFTri.kr((aL[i].sqrt - [aL[(i+1).mod(4)].sqrt * ~weights.kr(1,(i+1).mod(4)), aL[(i-1).mod(4)].sqrt * ~weights.kr(1,(i-1).mod(4))].mean))}).clip(-10.0, 10.0);
				var b = LPF.ar(pickIn * (Integrator.ar(sIn.abs,0.99 - freqDiff.clip(1,300).reciprocal.linlin(0,1.0, 0.01, 0.8).sqrt).reciprocal), 8000);
				var c = XFade2.ar(pickIn, b, a.linlin(0.0001,0.4, 1.0, -1.0));
				var d = XFade2.ar(sIn, c, A2K.kr(AnalogIn.ar(DC.ar(1))).round(0.001).linlin(0.0, 0.827, -1.0, 1.0));
				var comb = LPF.ar(4.collect({ |i|Mix(CombC.ar(d[i], 0.2, (RMS.kr(freqZero[i], 0.1) * [7/4, 9/8, 11/8, 13/8]).reciprocal,ampDiff.sqrt.linlin(0.0, 1, 0.5, -0.5)))} * 0.25), 4000);
				var out = LPF.ar(XFade2.ar(Mix(d), Mix(comb), AnalogIn.ar(2).round(0.001).linlin(0.0, 0.827, -1.0,1.0)), 6000);
				var locOut = Out.kr(~syncAmps, sync);
				var phaseImpulse = Impulse.kr(~distRange.kr);
				var impulse = Impulse.kr(10);
				var counter = PulseCount.kr(impulse, ~countReset.kr);
				var phasor = Gate.kr(Phasor.kr(~trig.kr, (~distRange.kr / sr)),~phasorGate.kr) ;
				SendReply.kr(impulse, '/phasor', a ++ [phasor, rmsQuick, counter]);
				Out.ar(0, Limiter.ar(LeakDC.ar(out)),lim);

		}).send(s)
	);

s.sync;

~exec = {
~trig.set(1); // Trigger phasor to reset
~phasorGate.set(1); // Open phasor gate
	"New Generation!".postln;
~distRange.set((d[\varSum][0] / d[\counter][0])); // Map distribution range to bus
	~variance.set((d[\variance])); // Map variance to bus

		4.collect({|i| d.add((\avgAmp_ ++ i).asSymbol -> (d[(\amp_ ++ i).asSymbol] / d[\counter][0]))}); //mean amp for each string

		4.collect({|i| d.add(\avgMainAmp -> (d[\mainAmp][0] / d[\counter][0]))}); //Mean overall amp

		4.collect({|i| d[\ampSync].put(i, (d[(\avgAmp_ ++ i).asSymbol][0] - (d[(\avgAmp_ ++ ((i-1).mod(4))).asSymbol][0] + d[(\avgAmp_++ ((i+1).mod(4))).asSymbol][0])/2).pow(2))}); // Mean variance for each string against its two neighbours for the past generation

		4.collect({|i| d[\ampSyncMono].put(i*2,(d[(\avgAmp_ ++ i).asSymbol][0] - d[(\avgAmp_ ++ ((i-1).mod(4))).asSymbol][0]).pow(2)); // Amp difference for string i against i - 1
			d[\ampSyncMono].put(i*2+1, (d[(\avgAmp_ ++ i).asSymbol][0] - d[(\avgAmp_ ++ (i+1).mod(4)).asSymbol][0]).pow(2)); // Amp difference for string i against i + 1
});

	d.add(\ampSyncMonoW -> (1 - d[\ampSyncMono]).normalizeSum); // Fitness landscape for individual syncs
	d.add(\ampSyncW -> (1 - d[\ampSync]).normalizeSum); // Fitness for overall sync
		8.collect({ |i| d[\fitnessLS].put(i, (d[\ampSyncMonoW][i] + d[\ampSyncW][((i.trunc(2))/2).asInt])/2)}); // Combined fitness landscape

4.collect({|i|
		d.add((\cut_ ++ i).asSymbol -> (rand(12) + 2)); //choose a random crossover point between 2 and 14
		~tempSS = ~searchSpace; // Make an array to store the temporary search space
			d.add((\choice_ ++(i * 2)).asSymbol -> ~searchSpace.wchoose(d[\fitnessLS])); // Choose the first parent index (weigthed probability based on fitness landscape)
		d[\tempSS].remove(d[(\choice_ ++ (i * 2)).asSymbol]); // Remove the first parent index from the temporary search space (so they don't mate with themselves)
			d.add((\choice_ ++ (i * 2 + 1)).asSymbol -> ~tempSS.choose);	// Choose the second parent index (weighted probability)
	}); // Creates an array of 8 values = 4 pairs of indexes for the parent weights

136.collect({|i|											// Mutation. Frequency = 1 / number of Genes
	if ((1/136).coin == true, { 							// If mutation happens
	var chromo = ((i.trunc(17)) / 17).asInt;
	var gene = i.mod(16);
	var biVal = d[(\parent_ ++ chromo).asSymbol][gene]; 	// Binary value of gene
	d[(\parent_ ++ chromo).asSymbol].put(gene, (biVal + 1).mod(2)); // Bit flip
	"Mutation!".postln;
})});


4.collect({|n|
	2.collect({|i| d.add((\parent_ ++ (i + (n * 2)) ++ "_0").asSymbol -> d[(\parent_ ++ d[(\choice_ ++ (i +  (n * 2))).asSymbol]).asSymbol].copyRange(0,d[(\cut_ ++ n).asSymbol]))}); // Archive the left part of the chromosome for each parent
	2.collect({|i| d.add((\parent_ ++ (i + (n * 2)) ++ "_1").asSymbol -> d[(\parent_ ++ d[(\choice_ ++ (i + (n * 2))).asSymbol]).asSymbol].copyRange(d[(\cut_ ++ n).asSymbol] + 1, 16))}) // Archive the right part of the chromosome for each parent
});

4.collect({|i|
			d.add((\parent_ ++ (i * 2)).asSymbol -> (d[(\parent_ ++ (i * 2) ++ "_0").asSymbol] ++ d[(\parent_ ++ (i * 2 + 1) ++ "_1").asSymbol])); // Generate new individual 1
	d.add((\parent_ ++ (i * 2 + 1)).asSymbol -> (d[(\parent_ ++ (i * 2 + 1) ++ "_0").asSymbol] ++ d[(\parent_ ++ (i * 2) ++ "_1").asSymbol]))});  // Generate new individual 2


	~weights.setn(8.collect({|i| (d[(\parent_++i).asSymbol].convertDigits(2)/100000 + 0.5) // Convert to binary digits. To do: use Grey code instead
	})); // Convert new generation of binary digits to weights between 0.5 and 1.5 and assign
		4.collect({|i| d[(\amp_ ++ i).asSymbol].put(0,0)}); // Reset individual amps
d[\mainAmp].put(0,0); //Reset main amp
d[\varSum].put(0,0); // Reset sum of variance
	~countReset.set(1);
	
				
				
					//"New gen born!".postln;
	};

	s.sync;
	~distRange.set(1);
	~weights.setn(8.collect({rrand(0.5, 1.5)}));
	x = Synth(\strings, target:s)
};



