
s = Server.default;

s.options.numAnalogInChannels = 8;
s.options.numAnalogOutChannels = 8;
s.options.numDigitalChannels = 16;

s.options.pgaGainLeft = 5;
s.options.pgaGainRight = 5;
s.options.dacLevel = 0;
s.options.adcLevel = 0;
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
~stringAmps = [2,2,0.8,1];
~syncAmps = Bus.control(s,4); // Amp for scaling incoming amplitudes
~syncFreq = Bus.control(s,4); // Frequency control for oscillators
~variance = Bus.control(s,1).set(0.0001); // Distrubtion range for sum of amplitudes
~weights = Bus.control(s,8).setn(1,1,1,1,1,1,1,1);
~phaseFreq = Bus.control(s,1).set(0.0001); // Frequency of phasor
~trig = Bus.control(s,1).set(1);
~distRange = Bus.control(s,1).set(0.0001); // Distrubtion range for sum of amplitudes
~countReset = Bus.control(s,1).set(0);
~phasorGate = Bus.control(s,1).set(1);
d = Dictionary.new; // Make a dictionary for holding all information about the GA
4.collect({|i| d.add((\amp_ ++ i).asSymbol -> [0])});
d.add(\mainAmp -> [0]);
d.add(\ampSync -> Array.fill(4,{0}));
d.add(\ampSyncMono -> Array.fill(8,{0}));
~searchSpace = 8.collect({|i| i}); // Make the searchspace (8 individual weights)
~tempSS = Array.fill(8, 0); // Temporary search space for storing parent indexes
8.collect({|i| d.add((\parent_ ++ i).asSymbol -> rand(100000).asBinaryDigits(17))}); // Generate 8 random weights of 17 binary digits each
d.add(\counter -> 1);
d.add(\variance -> 0);
d.add(\fitnessLS -> Array.fill(8,{rand(1.0)}).normalizeSum);


o = OSCFunc({|msg|
//	msg.postln;
//msg[7].postln;
var mean = msg[3..6].mean;
var distR = (4.collect({|i| (msg[i+3] - mean).pow(2)}).sum/4);
d[\mainAmp].put(0, d[\mainAmp][0]+msg[8]);
//d[\amp_0].postln;
4.collect({|i| d[(\amp_ ++ i).asSymbol].put(0,d[(\amp_ ++ i).asSymbol][0]+ msg[3+i])});
		if(msg[7] < 0.9,
		{~trig.set(0);
		~phasorGate.set(1); // Keep phasor gate open
			~countReset.set(0)},
		{
		~phasorGate.set(-1); // Keep phasor at 0.91
			if(distR > 0.8, // If distribution range is large
				{
					~distRange.set.(distR.linlin(0,3.0, 0.01,0.2)); // Map distribution range to bus
			//	~weights.setn(8.collect({rrand(0.2, 1.2)}).lag(0.2));
				~exec.value;
			})
		}
	)
},'/phasor');

s.waitForBoot {
	"Server Booted".postln;
	(
	SynthDef(\strings,
			{ |rmsLP = 4, filFrq =0.1, lim =0.9, oscMax =0.15, ampSc =33, ampDiffScale =1, dcScale = 1,  cD=0.2, cDecay= 0.2, ringTime =0.9 |
		var locIn = ~syncAmps.kr;
				var anIn = AnalogIn.ar(DC.ar(~anIns));
				var an2in = AnalogIn.ar(DC.ar(1)).linlin(0.002,0.827, 0.1,1.0).pow(2);
				var anScale = LinLin.ar(anIn, 0.002,0.827, 0.0,1.0).pow(2); // Scaled to 0-1
				var pickIn= SoundIn.ar(DC.ar([0,1,2,3]));
			var sIn = (pickIn*~stringAmps*locIn)*(anScale*ampSc); // Input scaled with controllers
				var freq = Pitch.kr(Mix(sIn), 440, 30)[0]; // Pitch detector
				var rmsFreq = RMS.kr(freq, 10);
				var rmsFreqSlow = RMS.kr(freq, 0.2);
				var freqDiff = RMS.kr(rmsFreq-rmsFreqSlow,1);
				var ab = RMS.kr(Mix(sIn), 0.2); // RMS amplitude follower, sum of all inputs
			var rmsQuick = RMS.kr(Mix(sIn), 10); // Fast RMS amplitude follower
				var ampDiff = RMS.kr((ab-rmsQuick)*ampDiffScale,4); // Difference between almost instantaneous amplitude and a running RMS
		//var a = RMS.kr(pickIn*an2in, 1); // RMS amplitude follower
		var a = RMS.kr(sIn,ampDiff.pow(2)); // RMS amplitude follower
	var sync = 4.collect({ |i| LFTri.kr(a[i]-[a[(i+1).mod(4)]*~weights.kr(1,(i+1).mod(4)), a[(i-1).mod(4)]*~weights.kr(1,(i-1).mod(4))].mean).range(0,1.0).clip(-10,10)}); 
				var clean = SoundIn.ar(DC.ar([0,1,2,3]),~stringAmps*anScale*ampSc); // Input scaled with controllers
				var	hpf = HPF.ar(XFade2.ar(Mix.new(clean), Mix.new(sIn),AnalogIn.ar(DC.ar(0)).linlin(0.002,0.827, -1.0,1.0)), 40);
				var comb = LPF.ar(CombC.ar(hpf, 0.2, `(rmsFreqSlow*[3/2,4/3, 7/4, 9/4, 11/4, 13/4, 17/8, 21/8]).reciprocal,rmsQuick.linlin(0.01,4.0, 0.01, 0.5).clip(0.01, 0.5)),5000) ;
				var comb2 = XFade2.ar(hpf, comb, (freqDiff.linlin(1,60.0, 2,0)-1));
				var out = XFade2.ar(hpf, comb2, AnalogIn.ar(DC.ar(2)).linlin(0.002,0.827, -1.0,1.0));
				var locOut = Out.kr(~syncAmps, sync);
				var phasor = Gate.kr(Phasor.kr(~trig.kr, ~distRange.kr),~phasorGate.kr) ;
				SendReply.kr(Impulse.kr(10), '/phasor', a ++ [phasor, rmsQuick]);
				Out.ar(1, Limiter.ar(LeakDC.ar(out)),lim, 0.001)
		}).send(s)
	);

s.sync;

~exec = {
~trig.set(1); // Trigger phasor to reset
~phasorGate.set(1); // Open phasor gate
	"New Generation!".postln;

	~variance.set((d[\variance].linlin(0,3.0, 0.01,0.2))); // Map variance to bus
			//	~weights.setn(8.collect({rrand(0.0, 1.0)}));

		4.collect({|i| d.add((\avgAmp_ ++ i).asSymbol -> (d[(\amp_ ++ i).asSymbol]/d[\counter]))}); //mean amp for each string

		4.collect({|i| d.add(\avgMainAmp -> (d[\mainAmp][0]/d[\counter]))}); //Mean overall amp

		4.collect({|i| d[\ampSync].put(i, (d[(\avgAmp_ ++ i).asSymbol][0] - (d[(\avgAmp_ ++ ((i-1).mod(4))).asSymbol][0] + d[(\avgAmp_++ ((i+1).mod(4))).asSymbol][0])/2).pow(2))}); //Variance for each string against its two neighbours

		4.collect({|i| d[\ampSyncMono].put(i*2,(d[(\avgAmp_ ++ i).asSymbol][0] - d[(\avgAmp_ ++ ((i-1).mod(4))).asSymbol][0]).pow(2)); // Amp difference for string i against i-1
			d[\ampSyncMono].put(i*2+1, (d[(\avgAmp_ ++ i).asSymbol][0] - d[(\avgAmp_ ++ (i+1).mod(4)).asSymbol][0]).pow(2));
// Amp difference for string i against i-1
});

	d.add(\ampSyncMonoW -> (1-d[\ampSyncMono]).normalizeSum); // Fitness landscape for individual syncs
	d.add(\ampSyncW -> (1-d[\ampSync]).normalizeSum); // Fitness for overall sync
		8.collect({ |i| d[\fitnessLS].put(i, (d[\ampSyncMonoW][i] + d[\ampSyncW][((i.trunc(2))/2).asInt])/2)}); // Combined fitness landscape

4.collect({|i|
		d.add((\cut_ ++ i).asSymbol -> (rand(12)+2)); //choose a random crossover point between 2 and 14
		~tempSS = ~searchSpace; // Make an array to store the temporary search space
			d.add((\choice_ ++(i * 2)).asSymbol -> ~searchSpace.wchoose(d[\fitnessLS])); // Choose the first parent index (weigthed probability based on fitness landscape)
		d[\tempSS].remove(d[(\choice_ ++ (i*2)).asSymbol]); // Remove the first parent index from the temporary search space (so they don't mate with themselves)
			d.add((\choice_ ++ (i*2+1)).asSymbol -> ~tempSS.choose);	// Choose the second parent index (weighted probability)
	}); // Creates an array of 8 values = 4 pairs of indexes for the parent weights

136.collect({|i| // Mutation. Frequency = 1/number of Genes
	if ((1/136).coin == true, { // If mutation happens
	var chromo = ((i.trunc(17))/17).asInt;
	var gene = i.mod(16);
	var biVal = d[(\parent_ ++ chromo).asSymbol][gene]; // Binary value of gene
	d[(\parent_ ++ chromo).asSymbol].put(gene, (biVal + 1).mod(2)); // Bit flip
	"Mutation!".postln;
})});


4.collect({|n|
	2.collect({|i| d.add((\parent_ ++ (i + (n*2)) ++ "_0").asSymbol -> d[(\parent_ ++ d[(\choice_ ++ (i +  (n*2))).asSymbol]).asSymbol].copyRange(0,d[(\cut_ ++ n).asSymbol]))}); // Archive the left part of the chromosome for each parent
	2.collect({|i| d.add((\parent_ ++ (i + (n*2)) ++ "_1").asSymbol -> d[(\parent_ ++ d[(\choice_ ++ (i + (n*2))).asSymbol]).asSymbol].copyRange(d[(\cut_ ++ n).asSymbol]+1, 16))}) // Archive the right part of the chromosome for each parent
});

4.collect({|i|
			d.add((\parent_ ++ (i*2)).asSymbol -> (d[(\parent_ ++ (i*2) ++ "_0").asSymbol] ++ d[(\parent_ ++ (i*2+1) ++ "_1").asSymbol])); // Generate new individual 1
	d.add((\parent_ ++ (i*2 + 1)).asSymbol -> (d[(\parent_ ++ (i*2+1) ++ "_0").asSymbol] ++ d[(\parent_ ++ (i*2) ++ "_1").asSymbol]))});  // Generate new individual 2


	~weights.setn(8.collect({|i| (d[(\parent_++i).asSymbol].convertDigits(2)/100000+0.5)
	})); // Convert new generation of binary digits to weights between 0.0 and 1 and assign
		4.collect({|i| d[(\amp_ ++ i).asSymbol].put(0,0)}); // Reset individual amps
d[\mainAmp].put(0,0); //Reset main amp
	~countReset.set(1);
	
				
				
					//"New gen born!".postln;
	};

	s.sync;
	~distRange.set(0.0001);
	~weights.setn(8.collect({rrand(0.5, 1.5)}));
	x = Synth(\strings, target:s)
};


