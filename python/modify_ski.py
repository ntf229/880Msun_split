import xml.etree.ElementTree as ET
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--filePath") # path to .ski file
parser.add_argument("--inc") # inclination angle (SKIRT parameter)
parser.add_argument("--numPhotons") # number of photon packages (SKIRT parameter)
parser.add_argument("--pixels") # number of pixels (square) for image (SKIRT parameter)
args = parser.parse_args()

tree = ET.parse(args.filePath)
root = tree.getroot()

d = {
	'FullInstrument/inclination' : str(args.inc)+'_deg',
    'FullInstrument/numPixelsX' : str(args.pixels),
    'FullInstrument/numPixelsY' : str(args.pixels),
	'MonteCarloSimulation/numPackets' : str(args.numPhotons),
}

for name, value in d.items():
	s = name.split('/')
	print('split:', s)
	print('length:', len(s))
	print(s[-1], value)

	for item in root.iter(s[0]):
		if len(s) == 2:
			item.set(s[-1], value.replace("_", " "))
		if len(s) == 3:
			for sub_item in item:
				sub_item.set(s[-1], value.replace("_", " "))

tree.write(args.filePath, encoding='UTF-8', xml_declaration=True)


