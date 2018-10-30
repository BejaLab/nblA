	function loadJSON(callback) {   

		var xobj = new XMLHttpRequest();
			xobj.overrideMimeType("application/json");
		xobj.open('GET', jsonFile, true);
		xobj.onreadystatechange = function () {
			  if (xobj.readyState == 4 && xobj.status == "200") {
				callback(xobj.responseText);
			  }
		};
		xobj.send(null);  
	}

	var consensus = {
		'%': { cutoffpercent: 0.60, residues: /[wlvimafcyhp]/i },
		'#': { cutoffpercent: 0.80, residues: /[wlvimafcyhp]/i },
		'-': { cutoffpercent: 0.50, residues: /[ed]/i },
		'+': { cutoffpercent: 0.60, residues: /[kr]/i },
		'g': { cutoffpercent: 0.50, residues: /[g]/i },
		'n': { cutoffpercent: 0.50, residues: /[n]/i },
		'q': { cutoffpercent: 0.50, residues: /[qe]/i },
		'p': { cutoffpercent: 0.50, residues: /[p]/i },
		't': { cutoffpercent: 0.50, residues: /[ts]/i },
		'A': { cutoffpercent: 0.85, residues: /[a]/i },
		'C': { cutoffpercent: 0.85, residues: /[c]/i },
		'D': { cutoffpercent: 0.85, residues: /[d]/i },
		'E': { cutoffpercent: 0.85, residues: /[e]/i },
		'F': { cutoffpercent: 0.85, residues: /[f]/i },
		'G': { cutoffpercent: 0.85, residues: /[g]/i },
		'H': { cutoffpercent: 0.85, residues: /[h]/i },
		'I': { cutoffpercent: 0.85, residues: /[i]/i },
		'K': { cutoffpercent: 0.85, residues: /[k]/i },
		'L': { cutoffpercent: 0.85, residues: /[l]/i },
		'M': { cutoffpercent: 0.85, residues: /[m]/i },
		'N': { cutoffpercent: 0.85, residues: /[n]/i },
		'P': { cutoffpercent: 0.85, residues: /[p]/i },
		'Q': { cutoffpercent: 0.85, residues: /[q]/i },
		'R': { cutoffpercent: 0.85, residues: /[r]/i },
		'S': { cutoffpercent: 0.85, residues: /[s]/i },
		'T': { cutoffpercent: 0.85, residues: /[t]/i },
		'V': { cutoffpercent: 0.85, residues: /[v]/i },
		'W': { cutoffpercent: 0.85, residues: /[w]/i },
		'Y': { cutoffpercent: 0.85, residues: /[y]/i },
	};

	var colors = {
		G: { color: '#F09048', conditions: "" },
		P: { color: '#C0C000', conditions: "" },
		T: { color: '#15C015', conditions: /[tST%#]/ },
		S: { color: '#15C015', conditions: /[tST#]/ },
		N: { color: '#15C015', conditions: /[nND]/ },
		Q: { color: '#15C015', conditions: /[qQE+KR]/ },
		W: { color: '#80A0F0', conditions: /[%#ACFHILMVWYPp]/ },
		L: { color: '#80A0F0', conditions: /[%#ACFHILMVWYPp]/ },
		V: { color: '#80A0F0', conditions: /[%#ACFHILMVWYPp]/ },
		I: { color: '#80A0F0', conditions: /[%#ACFHILMVWYPp]/ },
		M: { color: '#80A0F0', conditions: /[%#ACFHILMVWYPp]/ },
		A: { color: '#80A0F0', conditions: /[%#ACFHILMVWYPpTSsG]/ },
		F: { color: '#80A0F0', conditions: /[%#ACFHILMVWYPp]/ },
		C: { color: '#80A0F0', conditions: /[%#AFHILMVWYSPp]/ },
		H: { color: '#15A4A4', conditions: /[%#ACFHILMVWYPp]/ },
		Y: { color: '#15A4A4', conditions: /[%#ACFHILMVWYPp]/ },
		E: { color: '#C048C0', conditions: /[-DEqQ]/ },
		D: { color: '#C048C0', conditions: /[-DEnN]/ },
		K: { color: '#F01505', conditions: /[+KRQ]/ },
		R: { color: '#F01505', conditions: /[+KRQ]/ },
	};
	//var colors2 = {
	//	C: { color: '#FFC0CB', conditions: /[C]/ },
	//}

	window.addEventListener('load',function(){
		loadJSON(function(response) {
			var H = JSON.parse(response).height_arr;
			var picked = [];
			for (let i = 0; i < H.length; i++) {
				picked[i] = {};
				var matches = {};
				var cons = "";
				var sum = 0;
				for (let j = 0; j < H[i].length; j++) {
					var [ , aa, freq ] = H[i][j].match(/([A-Z]):([0-9.]+)/);
					freq = parseFloat(freq);
					sum += freq;
					Object.keys(consensus).forEach(function(key) {
						if (aa.match(consensus[key].residues)) {
							matches[key] = (matches[key] || 0) + freq;
						}
					});
				}
				Object.keys(matches).forEach(function(key) {
					if (matches[key]/sum >= consensus[key].cutoffpercent) cons = key;
				});
				Object.keys(colors).forEach(function(aa) {
					picked[i][aa] = (!colors[aa].conditions || cons.match(colors[aa].conditions)) ? colors[aa].color : "#C0C0C0";
				});
			}
			var texts = document.getElementsByTagName("text");
			var i = 0;
			var wasAA = false;
			for (let j = 0; j < texts.length; j++) {
				var aa = texts[j].innerHTML;
				if (aa.match(/^[A-Z]$/)) {
					texts[j].style.fill = picked[i][aa];
					wasAA = true;
				} else {
					if (wasAA) i++;
					wasAA = false;
				}
			};
		});
	});
