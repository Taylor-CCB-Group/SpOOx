from tifffile import TiffFile
import argsparse

parser = argparse.ArgumentParser(description='''Generate merged cellData.tab files so can run analyses at the source/sample level based on the files one level below
''', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--indir', dest='indir', default = "signalextraction",
                    help='initial directory level')
parser.add_argument('--excludeList', dest='excludeList', nargs='+',
                    help='list of dirs to exclude as a pattern if you have bad samples for example e.g. --excludeList AA_SAMPLE_11_ROI_1 BB_SAMPLE_4_ROI_12')
parser.add_argument('--infile', dest='infile', default = "cellData.tab",
                    help='the file it will look to merge')		
parser.add_argument('--outfile', dest='outfile', default = "mergecellData.tab",
                    help='the file to write out')			
parser.add_argument('--minarea', dest='minarea', default = "50", nargs='?', const=1, type=int,
                    help='minimum value for area of each cell')		
parser.add_argument('--maxarea', dest='maxarea', default = "300", nargs='?', const=1, type=int,
                    help='maximum value for area of each cell')		

args = parser.parse_args()

