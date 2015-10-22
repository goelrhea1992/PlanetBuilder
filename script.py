import sys
import os

input_path = sys.argv[1]
output_path = sys.argv[2]

files = os.listdir(input_path)

for conf in files:
	output_file = conf.split('.')[0] + "_output"
	command = 'java pb.sim.Simulator -g g9 --state ' + input_path + "/" + conf + ' > ' + output_path + output_file
	print command
	os.system(command)