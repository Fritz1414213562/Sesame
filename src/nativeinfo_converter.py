import sys
import numpy as np


def toml_force_field_local_bondlength(ninfo):

	retval = {"DNA": "", "Protein": ""}
	with open(ninfo, 'r') as fninfo:
		for line in fninfo:
			if line.startswith('bond'):
				tokens = line.split()
				isProtein = (len(tokens) == 13) and (tokens[12] == 'pp')
				isDNA     = len(tokens) == 12
				if not (isProtein or isDNA):
					print("Something wrong! The column number of 'bond' line is neither 12 or 13.\n{:s}".format(line), file=sys.stderr)
					sys.exit()
				params = {
					'record':         tokens[0],
					'iparam':         int(tokens[1]),
					'iunit1':         int(tokens[2]),
					'iunit2':         int(tokens[3]),
					'imp1':           int(tokens[4]),
					'imp2':           int(tokens[5]),
					'imp1un':         int(tokens[6]),
					'imp2un':         int(tokens[7]),
					'bd_nat':         float(tokens[8]),
					'factor_bd':      float(tokens[9]),
					'correct_bd_mgo': float(tokens[10]),
					'coef_bd':        float(tokens[11]),
				}
				if isProtein:
					iatom = params["imp1"] - 1
					jatom = params["imp2"] - 1
					v0    = params["bd_nat"]
					k     = params["coef_bd"]
					retval["Protein"] += f"{{indices = [{iatom:>4d}, {jatom:>4d}], v0 = {v0:10.6f}, k = {k:10.6f}}},\n"
				elif isDNA:
					iatom = params["imp1"] - 1
					jatom = params["imp2"] - 1
					v0    = params["bd_nat"]
					k     = params["coef_bd"]
					retval["DNA"] += f"{{indices = [{iatom:>4d}, {jatom:>4d}], v0 = {v0:10.6f}, k = {k:10.6f}}},\n"
	return retval


def toml_force_field_gocontact(ninfo):

	import constant
	retval = ""
	with open(ninfo, 'r') as fninfo:
		for line in fninfo:
			if line.startswith('contact'):
				tokens = line.split()
				isAICG = (len(tokens) == 13) and (tokens[12] == 'p-p')
				isGO   = len(tokens) == 12
				if not (isAICG or isGO):
					print("Something wrong! The column number of 'contact' line is neither 12 or 13.\n{:s}".format(line), file=sys.stderr)
					sys.exit()
				else:
					params = {
						'record':         tokens[0],
						'iparam':         int(tokens[1]),
						'iunit1':         int(tokens[2]),
						'iunit2':         int(tokens[3]),
						'imp1':           int(tokens[4]),
						'imp2':           int(tokens[5]),
						'imp1un':         int(tokens[6]),
						'imp2un':         int(tokens[7]),
						'go_nat':         float(tokens[8]),
						'factor_go':      float(tokens[9]),
						'dummy_mgo':      float(tokens[10]),
						'coef_go':        float(tokens[11]),
					}
					iatom = params["imp1"] - 1
					jatom = params["imp2"] - 1
					v0    = params["go_nat"]
					k     = params["coef_go"]
					if k > constant.ZERO_JUDGE:
						retval += f"{{indices = [{iatom:>4d}, {jatom:>4d}], v0 = {v0:10.6f}, k = {k:10.6f}}},\n"
	return retval


def toml_force_field_bondangle(ninfo, topol):
	import constant
	retval = {"DNA": "", "Protein": "", "Flexible": ""}
	with open(ninfo, 'r') as fninfo:
		for line in fninfo:
			if line.startswith('aicg13') or line.startswith('angl'):
				tokens = line.split()
				isProtein = len(tokens) == 16 and tokens[15] == 'ppp' and line.startswith('aicg13')
				isDNA     = len(tokens) == 14 and line.startswith('angl')
				if not (isProtein or isDNA):
				#	print("Something wrong! The column number of 'aicg13' line is neither 14 or 16.\n{:s}".format(line), file=sys.stderr)
				#	sys.exit()
					pass
				elif isProtein:
					params = {
						'record':             tokens[0],
						'iparam':             int(tokens[1]),
						'iunit1':             int(tokens[2]),
						'iunit2':             int(tokens[3]),
						'imp1':               int(tokens[4]),
						'imp2':               int(tokens[5]),
						'imp3':               int(tokens[6]),
						'imp1un':             int(tokens[7]),
						'imp2un':             int(tokens[8]),
						'imp3un':             int(tokens[9]),
						'aicg13_nat':         float(tokens[10]),
						'factor_aicg13':      float(tokens[11]),
						'correct_aicg13_mgo': float(tokens[12]),
						'coef_aicg13':        float(tokens[13]),
						'wid_aicg13_gauss':   float(tokens[14]),
						#'batype':             columns[15]
					}
					iatom =  params["imp1"] - 1
					jatom =  params["imp2"] - 1
					katom =  params["imp3"] - 1
					v0    =  params["aicg13_nat"]
					k     = -params["coef_aicg13"]
					sigma =  params["wid_aicg13_gauss"]
					res   =  topol.atom(jatom).residue.name
					if -k > constant.ZERO_JUDGE:
						retval["Protein"] += f'{{indices = [{iatom:>4d}, {katom:>4d}], v0 = {v0:10.6f}, k = {k:10.6f}, "σ" = {sigma:8.4f}}},\n'
						retval["Flexible"] += f'{{indices = [{iatom:>4d}, {jatom:>4d}, {katom:>4d}], k = 1.0, y = "y1_{res}", d2y = "y2_{res}"}},\n'
					else:
						retval["Flexible"] += f'{{indices = [{iatom:>4d}, {jatom:>4d}, {katom:>4d}], k = 1.0, y = "y1_{res}", d2y = "y2_{res}"}},\n'
				elif isDNA:
					params = {
						'record':             tokens[0],
						'iparam':             int(tokens[1]),
						'iunit1':             int(tokens[2]),
						'iunit2':             int(tokens[3]),
						'imp1':               int(tokens[4]),
						'imp2':               int(tokens[5]),
						'imp3':               int(tokens[6]),
						'imp1un':             int(tokens[7]),
						'imp2un':             int(tokens[8]),
						'imp3un':             int(tokens[9]),
						'ba_nat':             float(tokens[10]),
						'factor_ba':          float(tokens[11]),
						'correct_ba_mgo':     float(tokens[12]),
						'coef_ba':            float(tokens[13]),
					}
					iatom = params["imp1"] - 1
					jatom = params["imp2"] - 1
					katom = params["imp3"] - 1
					v0    = params["ba_nat"] * np.pi / 180.0
					k     = params["coef_ba"]
					retval["DNA"] += f'{{indices = [{iatom:>4d}, {jatom:>4d}, {katom:>4d}], k = {k:10.6f}, v0 = {v0:10.6f}}},\n'
	return retval


def toml_force_field_dihedralangle(ninfo, topol):
	import constant
	retval = {"3SPN2": "", "3SPN2C_1": "", "3SPN2C_2": "", "Protein": "", "Flexible": ""}
	with open(ninfo, 'r') as fninfo:
		for line in fninfo:
			if line.startswith('aicgdih') or line.startswith('dihd'):
				tokens     = line.split()
				isProtein  = len(tokens) == 18 and tokens[17] == 'pppp' and line.startswith('aicgdih')
				is3SPN2C_1 = len(tokens) == 18 and tokens[17] == 'N2P1' and line.startswith('dihd')
				is3SPN2C_2 = len(tokens) == 18 and tokens[17] == 'N2P2' and line.startswith('dihd')
				is3SPN2    = len(tokens) == 17 and line.startswith('dihd')
				if not (isProtein or is3SPN2C_1 or is3SPN2C_2 or is3SPN2):
				#	print("Something wrong! The column number of 'aicg13' line is not 18.\n{:s}".format(line), file=sys.stderr)
				#	sys.exit()
					pass
				elif isProtein:
					params = {
						'record':          tokens[0],
						'iparam':          int(tokens[1]),
						'iunit1':          int(tokens[2]),
						'iunit2':          int(tokens[3]),
						'imp1':            int(tokens[4]),
						'imp2':            int(tokens[5]),
						'imp3':            int(tokens[6]),
						'imp4':            int(tokens[7]),
						'imp1un':          int(tokens[8]),
						'imp2un':          int(tokens[9]),
						'imp3un':          int(tokens[10]),
						'imp4un':          int(tokens[11]),
						'dih_nat':         float(tokens[12]),
						'factor_aicg14':   float(tokens[13]),
						'correct_dih_mgo': float(tokens[14]),
						'coef_dih_gauss':  float(tokens[15]),
						'wid_dih_gauss':   float(tokens[16]),
						#'batype':             columns[15]
					}
					iatom =  params["imp1"] - 1
					jatom =  params["imp2"] - 1
					katom =  params["imp3"] - 1
					latom =  params["imp4"] - 1
					v0    =  params["dih_nat"] * np.pi / 180.0
					k     = -params["coef_dih_gauss"]
					sigma =  params["wid_dih_gauss"]
					jres  =  topol.atom(jatom).residue.name
					kres  =  topol.atom(katom).residue.name
					if -k > constant.ZERO_JUDGE:
						retval["Protein"] += f'{{indices = [{iatom:>4d}, {jatom:>4d}, {katom:>4d}, {latom:>4d}], v0 = {v0:10.6f}, k = {k:10.6f}, "σ" = {sigma:8.4f}}},\n'
						retval["Flexible"] += f'{{indices = [{iatom:>4d}, {jatom:>4d}, {katom:>4d}, {latom:>4d}], k = 1.0, coef = "{jres}_{kres}"}},\n'
					else:
						retval["Flexible"] += f'{{indices = [{iatom:>4d}, {jatom:>4d}, {katom:>4d}, {latom:>4d}], k = 1.0, coef = "{jres}_{kres}"}},\n'
				elif is3SPN2C_1:
					params = {
						'record':          tokens[0],
						'iparam':          int(tokens[1]),
						'iunit1':          int(tokens[2]),
						'iunit2':          int(tokens[3]),
						'imp1':            int(tokens[4]),
						'imp2':            int(tokens[5]),
						'imp3':            int(tokens[6]),
						'imp4':            int(tokens[7]),
						'imp1un':          int(tokens[8]),
						'imp2un':          int(tokens[9]),
						'imp3un':          int(tokens[10]),
						'imp4un':          int(tokens[11]),
						'dih_nat':         float(tokens[12]),
						'factor_dih':      float(tokens[13]),
						'correct_dih_mgo': float(tokens[14]),
						'coef_dih_1':      float(tokens[15]),
						'coef_dih_3':      float(tokens[16]),
						#'batype':             columns[15]
					}
					iatom = params["imp1"] - 1
					jatom = params["imp2"] - 1
					katom = params["imp3"] - 1
					latom = params["imp4"] - 1
					v0    = params["dih_nat"] * np.pi / 180.0
					if v0 < -np.pi:
						v0 += 2 * np.pi
					if v0 >  np.pi:
						v0 -= 2 * np.pi
					k     =  0.478011
					n     = 1
					retval["3SPN2C_1"] += f'{{indices = [{iatom:>4d}, {jatom:>4d}, {katom:>4d}, {latom:>4d}], v0 = {v0:>10.6f}, k = {k:>10.6f}, n = {n:>4d}}},\n'
				elif is3SPN2C_2:
					params = {
						'record':          tokens[0],
						'iparam':          int(tokens[1]),
						'iunit1':          int(tokens[2]),
						'iunit2':          int(tokens[3]),
						'imp1':            int(tokens[4]),
						'imp2':            int(tokens[5]),
						'imp3':            int(tokens[6]),
						'imp4':            int(tokens[7]),
						'imp1un':          int(tokens[8]),
						'imp2un':          int(tokens[9]),
						'imp3un':          int(tokens[10]),
						'imp4un':          int(tokens[11]),
						'dih_nat':         float(tokens[12]),
						'factor_dih':      float(tokens[13]),
						'correct_dih_mgo': float(tokens[14]),
						'coef_dih_1':      float(tokens[15]),
						'coef_dih_3':      float(tokens[16]),
						#'batype':             columns[15]
					}
					iatom = params["imp1"] - 1
					jatom = params["imp2"] - 1
					katom = params["imp3"] - 1
					latom = params["imp4"] - 1
					v0_g  = params["dih_nat"] * np.pi / 180.0 + np.pi
					v0_c  = v0_g - np.pi
					if v0_c < -np.pi:
						v0_c += 2 * np.pi
					if v0_c >  np.pi:
						v0_c -= 2 * np.pi
					k_c   =   0.478011
					k_g   =  -1.67304
					sigma =   0.3
					n     = 1
					retval["3SPN2C_2"] += f'{{indices = [{iatom:>4d}, {jatom:>4d}, {katom:>4d}, {latom:>4d}], Cosine = {{v0={v0_c:>10.6f}, k = {k_c:>10.6f}, n = {n:2d}}}, Gaussian = {{v0 = {v0_g:>10.6f}, k = {k_g:>10.6f}, sigma = {sigma:>10.6f} }} }},\n'
				elif is3SPN2:
					params = {
						'record':          tokens[0],
						'iparam':          int(tokens[1]),
						'iunit1':          int(tokens[2]),
						'iunit2':          int(tokens[3]),
						'imp1':            int(tokens[4]),
						'imp2':            int(tokens[5]),
						'imp3':            int(tokens[6]),
						'imp4':            int(tokens[7]),
						'imp1un':          int(tokens[8]),
						'imp2un':          int(tokens[9]),
						'imp3un':          int(tokens[10]),
						'imp4un':          int(tokens[11]),
						'dih_nat':         float(tokens[12]),
						'factor_dih':      float(tokens[13]),
						'correct_dih_mgo': float(tokens[14]),
						'coef_dih_1':      float(tokens[15]),
						'coef_dih_3':      float(tokens[16]),
						#'batype':             columns[15]
					}
					iatom = params["imp1"] - 1
					jatom = params["imp2"] - 1
					katom = params["imp3"] - 1
					latom = params["imp4"] - 1
					v0    = params["dih_nat"] * np.pi / 180.0
					k     =  -1.43403
					sigma =   0.3
					retval["3SPN2"] += f'{{indices = [{iatom:>4d}, {jatom:>4d}, {katom:>4d}, {latom:>4d}], v0 = {v0:>10.6f}, k = {k:>10.6f}, sigma = {sigma:>10.6f}  }},\n'
	return retval


def toml_force_field_local_base_params(topol):
	import constant
	retval = {"BaseStacking": "", "BasePair": "", "CrossStacking": ""}

	istrand = -1
	for ires in range(topol.n_residues):
		residue = topol.residue(ires)
		res     = residue.name
		group   = constant.GROUP_NAME[res]
		if group == 'DNA':
			base_name = res.strip('D')
			index_phos  = topol.select(f'resid {ires} and name == "DP"')
			index_sugar = topol.select(f'resid {ires} and name == "DS"')
			index_base  = topol.select(f'resid {ires} and name == "DB"')
			if len(index_phos) == 0:
				istrand += 1
				index_sugar = index_sugar[0]
				index_base  = index_base[0]
				retval["BaseStacking"] += f'{{strand = {istrand:>2d}, Base = "{base_name:1s}", S = {index_sugar:>4d}, B = {index_base:>4d}}},\n'
				retval["BasePair"] += f'{{strand = {istrand:>2d}, Base = "{base_name:1s}", S = {index_sugar:>4d}, B = {index_base:>4d}}},\n'
				retval["CrossStacking"] += f'{{strand = {istrand:>2d}, Base = "{base_name:1s}", S = {index_sugar:>4d}, B = {index_base:>4d}}},\n'
			else:
				index_phos  = index_phos[0]
				index_sugar = index_sugar[0]
				index_base  = index_base[0]
				retval["BaseStacking"] += f'{{strand = {istrand:>2d}, Base = "{base_name:1s}", P = {index_phos:>4d}, S = {index_sugar:>4d}, B = {index_base:>4d}}},\n'
				retval["BasePair"] += f'{{strand = {istrand:>2d}, Base = "{base_name:1s}", P = {index_phos:>4d}, S = {index_sugar:>4d}, B = {index_base:>4d}}},\n'
				retval["CrossStacking"] += f'{{strand = {istrand:>2d}, Base = "{base_name:1s}", P = {index_phos:>4d}, S = {index_sugar:>4d}, B = {index_base:>4d}}},\n'
	return retval


def toml_force_field_excluded_volume_params(topol, exv_scale):
	import constant

	retval = {"Total": "", "DNA": ""}
	for iatom in range(topol.n_atoms):
		atom  = topol.atom(iatom)
		name  = atom.name
		res   = atom.residue.name
		radius = ...
		if constant.GROUP_NAME[res] == "DNA":
			radius = constant.PARAM_EXV_SIGMA[res] if name == "DB" else constant.PARAM_EXV_SIGMA[name]
		else:
			radius = constant.PARAM_EXV_SIGMA[res]
		radius = 0.5 * exv_scale * radius
		retval["Total"] += f'{{index = {iatom:>4}, radius = {radius:>10.6f}}},\n'
		if constant.GROUP_NAME[res] == "DNA":
			atom_type = constant.ATOM_TYPE_DNA[res][name]
			retval["DNA"] += f'{{index = {iatom:>4d},  kind = "{atom_type}"}},\n'
	return retval


def toml_force_field_debye_huckel_electrostatics_params(topol, respac_charges = None):

	retval = {"InterChain":"", "IntraDNA":""}
	charges = {}
	intra_dnacharges = {}
	for iatom in range(topol.n_atoms):
		atom = topol.atom(iatom)
		name = atom.name
		res  = atom.residue.name
		if res in ("ARG", "LYS"):
			charges[iatom] = 1.0
		elif res in ("ASP", "GLU"):
			charges[iatom] = -1.0
		elif res in ("DA", "DT", "DG", "DC") and name == "DP":
			charges[iatom] = -1.0
			intra_dnacharges[iatom] = -0.6

	if respac_charges != None:
		for respac_idx in respac_charges:
			charges[respac_idx] = respac_charges[respac_idx]

	charges = sorted(charges.items(), key = lambda item: item[0])
	intra_dnacharges = sorted(intra_dnacharges.items(), key = lambda item: item[0])

	for charge in charges:
		idx = charge[0]
		q   = charge[1]
		retval["InterChain"] += f"{{index = {idx:>4d}, charge = {q:>10.4f}}},\n"
	for charge in intra_dnacharges:
		idx = charge[0]
		q   = charge[1]
		retval["IntraDNA"] += f"{{index = {idx:>4d}, charge = {q:>10.4f}}},\n"

	return retval


def toml_force_field_hydrogen_bonding_params(ninfo, topol, pdns_scale):
	import constant

	retval = ""
	for ires in range(topol.n_residues):
		residue = topol.residue(ires)
		res     = residue.name
		group   = constant.GROUP_NAME[res]
		if group == "DNA":
			index_phos  = topol.select(f'resid {ires} and name == "DP"')
			index_sugar = topol.select(f'resid {ires} and name == "DS"')
			index_base  = topol.select(f'resid {ires} and name == "DB"')
			if len(index_phos) == 0:
				continue
			else:
				index_phos  = index_phos[0]
				index_sugar = index_sugar[0]
				index_base  = index_base[0]
				retval += f'{{index = {index_phos:>4d}, S3 = {index_sugar:>4d}, kind = "DNA"}},\n'

	with open(ninfo, 'r') as fninfo:
		for line in fninfo:
			if line.startswith('pdns'):
				tokens = line.split()
				params = {
					'record':              tokens[0],
					'iparam':          int(tokens[1]),
					'ichain':          int(tokens[2]),
					'ires_global':     int(tokens[3]),
					'ires_local':      int(tokens[4]),
					'distance':        float(tokens[5]),
					'theta':           float(tokens[6]),
					'phi':             float(tokens[7]),
					'scale_factor':    float(tokens[8]),
				}
				index   = params['ires_global'] - 1
				index_n = topol.select(f'index {index-1} and name == "CA"')
				index_c = topol.select(f'index {index+1} and name == "CA"')
				if (len(index_n) != 0) and (len(index_c) != 0):
					index_n  = index_n[0]
					index_c  = index_c[0]
					distance = params['distance']
					theta    = params['theta']
					phi      = params['phi']
					scale_factor = -params['scale_factor'] * pdns_scale
					retval += f'{{index = {index:>4d}, kind = "Protein", PN = {index_n:>4d}, PC = {index_c:>4d}, k = {scale_factor:>4.2f}, r0 = {distance:>8.4f}, theta0 = {theta:>8.4f}, phi0 = {phi:>8.4f}}},\n'
	return retval
