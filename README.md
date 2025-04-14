# Sesame
Say 'Open Sesame' and you'll find a key for OpenCafeMol simulations in your hand.
## Usage
1. Download source codes.
2. Prepare a coarse-grained force-field parameter file (.ninfo) and a coarse-grained structure file (.pdb) using CafeMol. 
3. Generate a simulation system by running subcommand 'system'.
	Prepare a toml-format file which describes group-names. The group-names are used to make ignore list of forcefields for a pair potential.
	``` group_ids.toml:toml
	[group]
	indices = [
	{name = "Protein1", start =   0, end =  99},
	{name = "Protein2", start = 100, end = 300},
	{name = "Protein3", start = 301, end = 524},
	{name =      "DNA", start = 525, end = 900},
	....,
	]
	```
	``` terminal: python
	python3 <path-to-sesame>/src/main.py system --topol <path-to-coarse-grained-structure (.pdb)> --output <path-to-output (.toml)> --group <path-to-group-id>/group_ids.toml --temp [<temperature>, default = 300] --conc [<ionic-strength>, default = 0.3]
	```
4. Generate a set of forcefield parameters by running subcommand 'forcefield'.
	Prepare a toml-format file which describes pairs of group to ignore their electrostatic interactions.
	``` ignore_list.toml:toml
	[ignore]
	ele.inter = [["Protein1", "Protein2"], ["Protein2", "Protein3"], ....]
	ele.intra = ["Protein1", "Protein2", ....]
	ele_dna.inter = [["DNA", "Protein1"], ....]
	ele_dna.intra = ["Protein1", ....]
	```
	``` terminal: python
	python3 <path-to-sesame>/src/main.py forcefield --pdb <path-to-coarse-grained-structure (.pdb)> --ninfo <path-to-coarse-grained-forcefield-parameter (.ninfo)> --output <path-to-output (.toml)> --ignore [<path-to-ignore-list>/ignore_list.toml, default = None] --respac [<path-to-respac-parameter (.inp)>, default = None] --exv [float, default = 1.1] --pdns [float, default = 1.0] --nucltype ['2' or '2c', default = 2c] [--use_cap]
	```

