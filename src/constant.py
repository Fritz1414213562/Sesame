import math
GROUP_NAME = {
	"ALA": "PRO", "ARG": "PRO", "ASN": "PRO", "ASP": "PRO",
	"CYS": "PRO", "GLN": "PRO", "GLU": "PRO", "GLY": "PRO",
	"HIS": "PRO", "ILE": "PRO", "LEU": "PRO", "LYS": "PRO",
	"MET": "PRO", "PHE": "PRO", "PRO": "PRO", "SER": "PRO",
	"THR": "PRO", "TRP": "PRO", "TYR": "PRO", "VAL": "PRO",
	"OTH": "PRO",
	"DA" : "DNA", "DT" : "DNA", "DG" : "DNA", "DC" : "DNA",
	"RA" : "RNA", "RU" : "RNA", "RG" : "RNA", "RC" : "RNA",
}
PARAM_EXV_SIGMA = { # From cafemol exv.para
	# Protein
	"ALA": 3.943, "ARG": 4.166, "ASN": 4.328, "ASP": 4.419, "CYS": 4.449,
	"GLN": 3.955, "GLU": 4.316, "GLY": 3.448, "HIS": 4.696, "ILE": 4.352,
	"LEU": 4.342, "LYS": 4.53,  "MET": 4.424, "PHE": 4.403, "PRO": 4.046,
	"SER": 3.927, "THR": 4.036, "TRP": 4.476, "TYR": 4.135, "VAL": 4.44,
	"OTH": 4.00,

	# DNA
	"DP": 3.112, "DS": 3.689, "DA": 4.341, "DG": 4.631, "DC": 4.95,  "DT": 5.142,

	# RNA
	"P":  4.000, "S":  4.000, "Ab": 4.000, "Gb": 4.000, "Cb": 4.00 , "Ub": 4.000,
	#"P":  3.112, "S":  3.689, "Ab": 4.500, "Gb": 4.500, "Cb": 4.50 , "Ub": 4.500,

}
ATOM_TYPE_DNA = {
	"DA": {"DP": "P", "DS": "S", "DB": "A"},
	"DT": {"DP": "P", "DS": "S", "DB": "T"},
	"DG": {"DP": "P", "DS": "S", "DB": "G"},
	"DC": {"DP": "P", "DS": "S", "DB": "C"},
}
ATOM_MASS = {
	# Protein (from cafemol general.para)
	'ALA':	71.09, 'ARG': 156.19, 'ASN': 114.11, 'ASP': 115.09,
	'CYS': 103.15, 'GLN': 128.14, 'GLU': 129.12, 'GLY':  57.05,
	'HIS': 137.14, 'ILE': 113.16, 'LEU': 113.16, 'LYS': 128.17,
	'MET': 131.19, 'PHE': 147.18, 'PRO':  97.12, 'SER':  87.08,
	'THR': 101.11, 'TRP': 186.21, 'TYR': 163.18, 'VAL':  99.14,
	'OTH': 100.0,

	# DNA (from cafemol general.para, original and 3SPN2 original paper [2013 Hinckley JCP])
	'DP':  94.97, 'DS': 83.11, 'DA': 134.1,
	'DG': 150.1,  'DC': 110.1, 'DT': 125.1,

	# RNA (from cafemol general.para, original and Hori-san's paper [2012 Hori JCTC])
	'P':  62.97, 'S': 131.11, 'Ab': 134.1,
	'Gb': 150.1,  'Cb': 110.1, 'Ub': 111.1,

	# DNA (from Open3SPN2 paper, table12)
	# 'DP':  94.9696, 'DS': 83.1104,  'DA': 134.122, 
	# 'DG': 150.1214, 'DC': 110.0964, 'DT': 125.1078,
}
RNA_BOND_TYPE = {
	"PS", "SP", "SA", "SU", "SG", "SC", "SR", "SY", "SN"
}
RNA_ANGLE_TYPE = {
	"ASP", "USP", "GSP", "CSP", "RSP", "YSP", "NSP",
	"PSP", "SPS",
	"PSA", "PSU", "PSG", "PSC", "PSR", "PSY", "PSN",
}
RNA_DIHEDRAL_TYPE = {
	"ASPS", "USPS", "GSPS", "CSPS", "RSPS", "YSPS", "NSPS",
	"PSPS", "SPSP",
	"SPSA", "SPSU", "SPSG", "SPSC", "SPSR", "SPSY", "SPSN",
}
NM2ANGSTROM = 10.0
KJPERKCAL = 4.1868
CAFETIME = math.sqrt(1.0 / KJPERKCAL) * 0.1 # ps/cafe-time
ZERO_JUDGE = 1e-6
