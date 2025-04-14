

def read_respac(fname):

	retval = dict()
	with open(fname, 'r') as fin:
		for line in fin:
			if line.startswith('CHARGE_CHANGE'):
				tokens = line.split()
				index = int(tokens[1]) - 1
				q     = float(tokens[2])
				retval[index] = q
	return retval
