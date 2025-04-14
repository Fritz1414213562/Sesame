

def run(args):
	main(args.topol, args.traj, args.velo, args.output, args.group, args.rst, args.temp, args.conc)


def main(pdb, traj, velo, output, group, rst_frame, temperature, concentration):
	import mdtraj
	import constant
	import toml

	fp = open(output, 'w')
	fp.write("[[systems]]\n")
	fp.write("attributes.temperature = {:5.2f}\n".format(temperature))
	fp.write("attributes.ionic_strength = {:5.2f}\n".format(concentration))

	# Particles
	fp.write("particles = [\n")

	trajectory = mdtraj.load_dcd(traj, top = pdb)
	topology = trajectory.topology
	velocity = mdtraj.load_dcd(velo, top = pdb)
	meta = toml.load(open(group))
	group_indices = meta["group"]["indices"]

	for iatom in range(trajectory.n_atoms):
		atom = topology.atom(iatom)
		name = atom.name
		res  = atom.residue.name
		x	 = trajectory.xyz[rst_frame, iatom, 0] * constant.NM2ANGSTROM
		y	 = trajectory.xyz[rst_frame, iatom, 1] * constant.NM2ANGSTROM
		z	 = trajectory.xyz[rst_frame, iatom, 2] * constant.NM2ANGSTROM
		vx   = velocity.xyz[rst_frame, iatom, 0] * constant.NM2ANGSTROM
		vy   = velocity.xyz[rst_frame, iatom, 1] * constant.NM2ANGSTROM
		vz   = velocity.xyz[rst_frame, iatom, 2] * constant.NM2ANGSTROM
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


