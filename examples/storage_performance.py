# Example Storage Performance
# Michael Wollensack METAS - 15.03.2019

print('Example Storage Performance')
print('Begin')

import os

from metas_unclib import *
use_linprop()

n = 20
fp_bin = 'testdata.bin'
fp_xml = 'testdata.xml'

# Create Random Data
data0 = ufloatarray(np.random.rand(n, n), np.eye(n * n))
data = ulinalg.inv(data0)

# Binary Serialization
print('Binaray serialization - Save file')
ustorage.save_binary_file(data, fp_bin)
print('Binaray serialization - Load file')
bin_data = ustorage.load_binary_file(fp_bin)
os.remove(fp_bin);

# Xml Serialization
print('Xml serialization - Save file')
ustorage.save_xml_file(data, fp_xml)
print('Xml serialization - Load file')
xml_data = ustorage.load_xml_file(fp_xml)
os.remove(fp_xml);

print('End')
