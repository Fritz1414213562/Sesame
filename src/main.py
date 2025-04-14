

def main():
	import argparse
	import make_forcefield as ff
	import make_system     as syst
	import restart_system  as rest

	parser = argparse.ArgumentParser(description = "Preparation of input files for OpenAICG2+ based on CafeMol input files")
	subparsers = parser.add_subparsers()

## ForceField
	parser_ff = subparsers.add_parser('forcefield', help = 'see `forcefield -h`')
	parser_ff.add_argument("--ninfo", '-n', nargs = '*', required = True)
	parser_ff.add_argument("--pdb", '-p', required = True)
	parser_ff.add_argument("--output", '-o', required = True)
	parser_ff.add_argument("--respac", '-q')
	parser_ff.add_argument("--ignore", '-i')
	parser_ff.add_argument("--exv", '-x', type = float, default = 1.1)
	parser_ff.add_argument("--pdns", '-y', type = float, default = 1.0)
	parser_ff.add_argument("--nucltype", '-t', default = "2c")
	parser_ff.add_argument("--use_cap", '-c', action = "store_true")
	parser_ff.set_defaults(handler = ff.run)

## System

	parser_sys = subparsers.add_parser('system', help = 'see `system -h`')
	parser_sys.add_argument("--topol", '-p', required = True)
	parser_sys.add_argument("--output", '-o', required = True)
	parser_sys.add_argument("--group", '-g', required = True)
	parser_sys.add_argument("--velo", '-v')
	parser_sys.add_argument("--rst", '-r', type = int, default = -1)
	parser_sys.add_argument("--temp", '-t', type = float, default = 300.0)
	parser_sys.add_argument("--conc", '-c', type = float, default = 0.3)
	parser_sys.set_defaults(handler = syst.run)

## Restart

	parser_res = subparsers.add_parser('restart', help = 'see `system -h`')
	parser_res.add_argument("--topol", '-p', required = True)
	parser_res.add_argument("--traj", '-j', required = True)
	parser_res.add_argument("--velo", '-v', required = True)
	parser_res.add_argument("--output", '-o', required = True)
	parser_res.add_argument("--group", '-g', required = True)
	parser_res.add_argument("--rst", '-r', type = int, default = -1)
	parser_res.add_argument("--temp", '-t', type = float, default = 300.0)
	parser_res.add_argument("--conc", '-c', type = float, default = 0.3)
	parser_res.set_defaults(handler = rest.run)

## Parse command-line arguments and run

	args = parser.parse_args()
	if hasattr(args, 'handler'):
		args.handler(args)
	else:
		parser.print_help()

if __name__ == "__main__":
	main()
