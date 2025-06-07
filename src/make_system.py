

def run(args):
	main(args.pdb, args.output, args.group, args.velo, args.rst, args.temp, args.conc)


def main(pdb, output, group, velo, rst_frame, temperature, concentration):
	import mdtraj
	import constant
	import toml

	fp = open(output, 'w')
	fp.write("[[systems]]\n")
	fp.write("attributes.temperature = {:5.2f}\n".format(temperature))
	fp.write("attributes.ionic_strength = {:5.2f}\n".format(concentration))

	# Particles
	fp.write("particles = [\n")

	traj = mdtraj.load_pdb(pdb)
	topo = traj.topology
	isVelocityLoaded = not (velo == None)
	velocity = mdtraj.load_dcd(velo, top = pdb) if isVelocityLoaded else None
	meta = toml.load(open(group))
	group_indices = meta["group"]["indices"]
#	vxyz = None if not velocity else velocity.xyz[rst_frame, :, :] * constant.NM2ANGSTROM * constant.CAFETIME
	vxyz = velocity.xyz[rst_frame, :, :] * constant.NM2ANGSTROM if isVelocityLoaded else None

	for iatom in range(traj.n_atoms):
		atom = topo.atom(iatom)
		name = atom.name
		res  = atom.residue.name
		x	 = traj.xyz[0, iatom, 0] * constant.NM2ANGSTROM
		y	 = traj.xyz[0, iatom, 1] * constant.NM2ANGSTROM
		z	 = traj.xyz[0, iatom, 2] * constant.NM2ANGSTROM
		vx   = 0.0 if not isVelocityLoaded else vxyz[iatom, 0]
		vy   = 0.0 if not isVelocityLoaded else vxyz[iatom, 1]
		vz   = 0.0 if not isVelocityLoaded else vxyz[iatom, 2]
		group = ...
		for items in group_indices:
			if items["start"] <= iatom <= items["end"]:
				group = items["name"]

		atom_type = ...
		mass = ...
		if group == 'DNA':
			atom_type = constant.ATOM_TYPE_DNA[res][name];
			mass = constant.ATOM_MASS[res] if name == "DB" else constant.ATOM_MASS[name]
		else:
			atom_type = name;
			mass = constant.ATOM_MASS[res]
		fp.write(f'{{m={mass:>6.2f}, pos=[{x:>8.3f}, {y:>8.3f}, {z:>8.3f}], vel=[{vx:>8.5f}, {vy:>8.5f}, {vz:>8.5f}], name="{atom_type}", group="{group}"}},\n')

	fp.write("]\n")
	fp.write('\n')
	fp.close()


if __name__ == "__main__":
	import argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("--pdb", '-p', required = True)
	parser.add_argument("--group", '-g', required = True)
	parser.add_argument("--output", '-o', required = True)
	parser.add_argument("--velo", '-v')
	parser.add_argument("--rst", '-f', type = int, default = -1)
	parser.add_argument("--temp", '-t', type = float, default = 300.0)
	parser.add_argument("--conc", '-c', type = float, default = 0.3)
	args = parser.parse_args()
	main(args.pdb, args.output, args.group, args.velo, args.rst, args.temp, args.conc)
