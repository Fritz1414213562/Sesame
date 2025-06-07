from nativeinfo_converter import *
from respac_reader import *
import numpy as np
import mdtraj
import toml
import sys


def run(args):
	main(args.ninfo, args.pdb, args.respac, args.exv, args.pdns, args.ignore, args.nucltype, args.use_cap, args.adjust_exclusion, args.output)


def main(ninfos, pdb, respac, exv_scale, pdns_scale, ignore, nucltype, use_capped_go, adjust_exclusion, out):

	bond_params = ...
	go_params = ...
	angl_params = ...
	dihd_params = ...
	pdns_params = ...
	isFirstNinfo = True
	nucleotide_types = ["2c", "2"]
	nucl_ff_name = ...

	if not (nucltype in nucleotide_types):
		print("The nucleotide type {:s} is not supported.\nSupported types are '2' and '2c'".format(nucltype), file = sys.stderr)
		sys.exit()
	elif nucltype == "2":
		nucl_ff_name = "3SPN2"
	elif nucltype == "2c":
		nucl_ff_name = "3SPN2C"
	else:
		print("The forcefield name for nucleotide type {:s} is not undefined.".format(nucltype), file = sys.stderr)
		sys.exit()

	topol = mdtraj.load_pdb(pdb).topology
	respac_charges = read_respac(respac) if respac != None else None

	base_params = toml_force_field_local_base_params(topol)
	exv_params  = toml_force_field_excluded_volume_params(topol, exv_scale)
	ele_params  = toml_force_field_debye_huckel_electrostatics_params(topol, respac_charges)

	for ninfo in ninfos:
		if isFirstNinfo:
			bond_params = toml_force_field_local_bondlength(ninfo)
			go_params   = toml_force_field_gocontact(ninfo)
			angl_params = toml_force_field_bondangle(ninfo, topol)
			dihd_params = toml_force_field_dihedralangle(ninfo, topol)
			pdns_params = toml_force_field_hydrogen_bonding_params(ninfo, topol, pdns_scale)
			isFirstNinfo = False
		else:
			tmp_bond_params = toml_force_field_local_bondlength(ninfo)
			tmp_go_params   = toml_force_field_gocontact(ninfo)
			tmp_angl_params = toml_force_field_bondangle(ninfo, topol)
			tmp_dihd_params = toml_force_field_dihedralangle(ninfo, topol)
			tmp_pdns_params = toml_force_field_hydrogen_bonding_params(ninfo, topol, pdns_scale)

			bond_params["DNA"]      += tmp_bond_params["DNA"]
			bond_params["Protein"]  += tmp_bond_params["Protein"]
			bond_params["RNA"]      += tmp_bond_params["RNA"]
			go_params               += tmp_go_params
			angl_params["DNA"]      += tmp_angl_params["DNA"]
			angl_params["Protein"]  += tmp_angl_params["Protein"]
			angl_params["Flexible"] += tmp_angl_params["Flexible"]
			angl_params["RNA"]      += tmp_angl_params["RNA"]
			dihd_params["3SPN2C_1"] += tmp_dihd_params["3SPN2C_1"]
			dihd_params["3SPN2C_2"] += tmp_dihd_params["3SPN2C_2"]
			dihd_params["3SPN2"]    += tmp_dihd_params["3SPN2"]
			dihd_params["Protein"]  += tmp_dihd_params["Protein"]
			dihd_params["Flexible"] += tmp_dihd_params["Flexible"]
			dihd_params["RNA_1"]    += tmp_dihd_params["RNA_1"]
			dihd_params["RNA_2"]    += tmp_dihd_params["RNA_2"]
			pdns_params             += tmp_pdns_params
	
	with open(out, 'w') as fout:
		## Local
		# Bond length
		fout.write('[[forcefields]]\n')
		if len(bond_params["Protein"]) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "BondLength"\n')
			fout.write('potential   = "Harmonic"\n')
			fout.write('topology    = "bond"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(bond_params["Protein"]))
		if len(bond_params["DNA"]) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "BondLength"\n')
			fout.write('potential   = "3SPN2Bond"\n')
			fout.write('topology    = "bond"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(bond_params["DNA"]))
		if len(bond_params["RNA"]) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "BondLength"\n')
			fout.write('potential   = "Harmonic"\n')
			fout.write('topology    = "bond"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(bond_params["RNA"]))
		# Go-contact
		if len(go_params) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "BondLength"\n')
			fout.write('potential   = "GoContact"\n') if not use_capped_go else fout.write('potential   = "CappedGoContact"\n')
			fout.write('topology    = "contact"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(go_params))
		# Angle
		if len(angl_params["Protein"]) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "BondLength"\n')
			fout.write('potential   = "Gaussian"\n')
			fout.write('topology    = "none"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(angl_params["Protein"]))
		if len(angl_params["Flexible"]) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "BondAngle"\n')
			fout.write('potential   = "FlexibleLocalAngle"\n')
			fout.write('topology    = "none"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(angl_params["Flexible"]))
		if len(angl_params["DNA"]) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "BondAngle"\n')
			fout.write('potential   = "Harmonic"\n')
			fout.write('topology    = "angle"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(angl_params["DNA"]))
		if len(angl_params["RNA"]) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "BondAngle"\n')
			fout.write('potential   = "Harmonic"\n')
			fout.write('topology    = "angle"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(angl_params["RNA"]))
		# Dihedral angle
		if len(dihd_params["Protein"]) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "DihedralAngle"\n')
			fout.write('potential   = "Gaussian"\n')
			fout.write('topology    = "none"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(dihd_params["Protein"]))
		if len(dihd_params["Flexible"]) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "DihedralAngle"\n')
			fout.write('potential   = "FlexibleLocalDihedral"\n')
			fout.write('topology    = "none"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(dihd_params["Flexible"]))
		if nucl_ff_name == "3SPN2C":
			if len(dihd_params["3SPN2C_1"]) < 1:
				print("Warning: No dihedral parameters for cosine-shaped potentials in 3SPN.2C", file = sys.stderr)
			else:
				fout.write('[[forcefields.local]]\n')
				fout.write('interaction = "DihedralAngle"\n')
				fout.write('potential   = "Cosine"\n')
				fout.write('topology    = "dihedral"\n')
				fout.write('parameters = [\n{:s}]\n\n'.format(dihd_params["3SPN2C_1"]))
			if len(dihd_params["3SPN2C_2"]) < 1:
				print("Warning: No dihedral parameters for cosine- and gaussian-shaped potentials in 3SPN.2C", file = sys.stderr)
			else:
				fout.write('[[forcefields.local]]\n')
				fout.write('interaction = "DihedralAngle"\n')
				fout.write('potential   = "Gaussian+Cosine"\n')
				fout.write('topology    = "dihedral"\n')
				fout.write('parameters = [\n{:s}]\n\n'.format(dihd_params["3SPN2C_2"]))
		elif nucl_ff_name == "3SPN2":
			if len(dihd_params["3SPN2"]) < 1:
				print("Warning: No dihedral parameters in 3SPN.2", file = sys.stderr)
			else:
				fout.write('[[forcefields.local]]\n')
				fout.write('interaction = "DihedralAngle"\n')
				fout.write('potential   = "Gaussian"\n')
				fout.write('topology    = "dihedral"\n')
				fout.write('parameters = [\n{:s}]\n\n'.format(dihd_params["3SPN2"]))
		if len(dihd_params["RNA_1"]) < 1:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "DihedralAngle"\n')
			fout.write('potential   = "Cosine"\n')
			fout.write('topology    = "dihedral"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(dihd_params["RNA_1"]))
		if len(dihd_params["RNA_2"]) < 1:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "DihedralAngle"\n')
			fout.write('potential   = "Cosine"\n')
			fout.write('topology    = "none"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(dihd_params["RNA_2"]))
		# Base
		if len(base_params["BaseStacking"]) > 0:
			fout.write('[[forcefields.local]]\n')
			fout.write('interaction = "3SPN2BaseStacking"\n')
			fout.write('potential   = "{:s}"\n'.format(nucl_ff_name))
			fout.write('topology    = "nucleotide"\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(base_params["BaseStacking"]))
		## Global
		# Base
		if len(base_params["BasePair"]) > 0:
			fout.write('[[forcefields.global]]\n')
			fout.write('interaction = "3SPN2BasePair"\n')
			fout.write('potential   = "{:s}"\n'.format(nucl_ff_name))
			if adjust_exclusion:
				fout.write('ignore.particles_within.nucleotide = 2\n')
			else:
				fout.write('ignore.particles_within.nucleotide = 3\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(base_params["BasePair"]))
		if len(base_params["CrossStacking"]) > 0:
			fout.write('[[forcefields.global]]\n')
			fout.write('interaction = "3SPN2CrossStacking"\n')
			fout.write('potential   = "{:s}"\n'.format(nucl_ff_name))
			if adjust_exclusion:
				fout.write('ignore.particles_within.nucleotide = 2\n')
			else:
				fout.write('ignore.particles_within.nucleotide = 3\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(base_params["CrossStacking"]))
		# Excluded Volume Interaction (w/o intra DNA)
		if len(exv_params["Total"]) > 0:
			fout.write('[[forcefields.global]]\n')
			fout.write('interaction = "Pair"\n')
			fout.write('potential = "ExcludedVolume"\n')
			fout.write('ignore.molecule = "Nothing"\n')
			fout.write('ignore.particles_within.bond = 3\n')
			fout.write('ignore.particles_within.contact = 1\n')
			fout.write('ignore.group.intra = ["DNA"]\n')
			fout.write('epsilon = 0.6\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(exv_params["Total"]))
		# Excluded Volume Interaction (for intra DNA)
		if len(exv_params["DNA"]) > 0:
			fout.write('[[forcefields.global]]\n')
			fout.write('interaction = "Pair"\n')
			fout.write('potential = "3SPN2ExcludedVolume"\n')
			fout.write('ignore.particles_within.angle      = 1\n')
			fout.write('ignore.particles_within.nucleotide = 1\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(exv_params["DNA"]))
		# Debye-Huckel Electrostatics for default
		if len(ele_params["InterChain"]) > 0:
			fout.write('[[forcefields.global]]\n')
			fout.write('interaction = "Pair"\n')
			fout.write('potential   = "DebyeHuckel"\n')
			fout.write('ignore.particles_within.bond = 3\n')
			fout.write('ignore.particles_within.nucleotide = 1\n')
			meta = toml.load(open(ignore)) if ignore else None
			if meta and ("ele" in set(meta["ignore"].keys())):
				if "intra" in set(meta["ignore"]["ele"].keys()):
					ignore_list = meta["ignore"]["ele"]["intra"]
					ignore_intra = '"' + '","'.join(ignore_list) + '"'
					fout.write("ignore.group.intra = [{:s}]\n".format(ignore_intra))
				if "inter" in set(meta["ignore"]["ele"].keys()):
					ignore_list = meta["ignore"]["ele"]["inter"]
					ignore_inter = ', '.join(['["{:s}", "{:s}"]'.format(pair[0], pair[1]) for pair in ignore_list])
					fout.write("ignore.group.inter = [{:s}]\n".format(ignore_inter))
			fout.write("parameters = [\n{:s}]\n\n".format(ele_params["InterChain"]))
		# Debye-Huckel Electrostatics for DNA-protein interaction
		if len(ele_params["IntraDNA"]) > 0:
			fout.write('[[forcefields.global]]\n')
			fout.write('interaction = "Pair"\n')
			fout.write('potential   = "DebyeHuckel"\n')
			fout.write('ignore.particles_within.nucleotide = 1\n')
			if meta and ("ele_dna" in set(meta["ignore"].keys())):
				meta = toml.load(open(ignore))
				if "intra" in set(meta["ignore"]["ele_dna"].keys()):
					ignore_list = meta["ignore"]["ele_dna"]["intra"]
					ignore_intra = '"' + '","'.join(ignore_list) + '"'
					fout.write("ignore.group.intra = [{:s}]\n".format(ignore_intra))
				if "inter" in set(meta["ignore"]["ele_dna"].keys()):
					ignore_list = meta["ignore"]["ele_dna"]["inter"]
					ignore_inter = ', '.join(['["{:s}", "{:s}"]'.format(pair[0], pair[1]) for pair in ignore_list])
					fout.write("ignore.group.inter = [{:s}]\n".format(ignore_inter))
			fout.write("parameters = [\n{:s}]\n\n".format(ele_params["IntraDNA"]))
		# Hydrogen-bonding
		if len(pdns_params) > 0:
			fout.write('[[forcefields.global]]\n')
			fout.write('interaction = "PDNS"\n')
			fout.write('potential   = "PDNS"\n')
			fout.write('sigma       = 1.0   \n')
			fout.write('delta       = 0.17453\n')
			fout.write('cutoff      = 10\n')
			fout.write('energy_unit = 0.593\n')
			fout.write('parameters = [\n{:s}]\n\n'.format(pdns_params))


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("--ninfo", '-n', nargs = '*', required = True)
	parser.add_argument("--pdb", '-p', required = True)
	parser.add_argument("--output", '-o', required = True)
	parser.add_argument("--respac", '-q')
	parser.add_argument("--ignore", '-i')
	parser.add_argument("--exv", '-x', type = float, default = 1.1)
	parser.add_argument("--pdns", '-y', type = float, default = 1.0)
	parser.add_argument("--nucltype", '-t', default = "2c")
	parser.add_argument("--use_cap", '-c', action = "store_true")
	parser.add_argument("--adjust_exclusion", '-a', action = "store_true")
	args = parser.parse_args()
	main(args.ninfo, args.pdb, args.respac, args.exv, args.pdns, args.ignore, args.nucltype, args.use_cap, args.adjust_exclusion, args.output)
