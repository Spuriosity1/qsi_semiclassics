import sys
import time

# usage: binary_file target.hhtemplate target.hh
assert len(sys.argv) > 3

with open(sys.argv[1],'rb') as f:
    data = []
    s = f.read(4)
    while len(s) == 4:
        data.append(int.from_bytes(s, byteorder='little'))
        s = f.read(4)

with open(sys.argv[2],'r') as ft:
    with open(sys.argv[3],'w') as fd:
        fd.write('''
// THIS FILE WAS AUTOMATICALLY GENERATED ON %s
// DO NOT EDIT IT MANUALLY, EDIT seed.hhtemplate AND REMAKE
''' % time.ctime())
        for line in ft:
            if line.startswith('/* __RANDOM_NUMBERS__ */'):
                for i, num in enumerate(data):
                    if i == 0:
                        fd.write('0x%x' % num)
                    else:
                        fd.write(',\n0x%x' % num)
                fd.write('\n')
            elif line.startswith('const size_t NUM_SEEDS ='):
                fd.write('const size_t NUM_SEEDS = %d;'% len(data))
            else:
                fd.write(line)


