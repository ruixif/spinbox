import numpy as np
from preutils import munck_parser, calibration_base, read_calibration, parse_arguments, folddata, guo_parser
from matplotlib import pyplot as plt
import csv

from os import listdir, mkdir
from os.path import isfile, join, isdir


class Spectrum:
    def __init__(self, source_file, file_type, calibration, **KWARGS):
        if file_type == "munck_cmu":
            self._filename, self._date, self._liner, self._temperature, self._field, self._veloscale, \
                 self._rawdataseq, self._wholeblock = munck_parser(source_file)
        elif file_type == "guo_cmu":
            self._filename, self._date, self._liner, self._temperature, self._field, self._veloscale, \
                 self._rawdataseq, self._wholeblock = guo_parser(source_file)
        else:
            print(file_type + "not available")
            
        self.calibration = calibration
        
    def calibrate(self):
        slope, ZVC, foldpoint = read_calibration(self.calibration, self._date, self._veloscale)
        print(slope, ZVC, foldpoint)
        if ('ff' in self._filename) or ('FF' in self._filename):
            #folded data
            self._count = self._rawdataseq
            channels = np.linspace(1, 256, num=256)
            velocity = [(each-ZVC) * slope for each in channels]
            self._velocity = velocity
        elif ('rr' in self._filename) or ('RR' in self._filename) or ('r' in self._filename):
            #unfolded data 
            channels = np.linspace(1, 512, num=512)
            aved_channel, self._count = folddata(foldpoint, channels, self._rawdataseq)
            self._velocity = [(each-ZVC) * slope for each in aved_channel]
        else:
            
            raise RuntimeError('Unrecognized file')
            
        
        self._errorbar = 100 / (np.sqrt(self._count[0]))
        baseline = sum(self._count[0:3])/3
        prenormal = [each - baseline for each in self._count]
        total_count_area = sum(prenormal) * slope * 10
        self._normal = [100 * each/total_count_area for each in prenormal]
        
    def plot(self, outpath, showflag = False):
        this_fig = plt.figure()
        plt.gca().invert_yaxis()
        this_fig.suptitle(self._wholeblock)
        plt.ylabel('absorption(%)')
        plt.xlabel('velocity(mm/s)')
        plt.errorbar(self._velocity, self._normal, self._errorbar)
        if showflag:
            this_fig.show()
        else:
            this_fig.savefig(outpath + '/' + self._filename + '.png')
        plt.close()
        
    def export(self, outpath=None):
        if outpath is None:
            outfile = self._filename + '.csv'
        else:
            outfile = outpath + '/' + self._filename + '.csv'
        with open(outfile, 'w') as out_ptr:
            out_ptr.write(self._wholeblock)
            writer = csv.writer(out_ptr)
            writer.writerows(zip(self._velocity, self._normal))

def main(args):
    args = parse_arguments()
    sourcepath = args.sourcepath
    plotflag = args.plotflag
    processedpath = args.processedpath
    onlyfiles = [join(sourcepath, files) for files in listdir(sourcepath) if isfile(join(sourcepath, files))]
    if isdir(processedpath):
        pass
    else:
        mkdir(processedpath)
    
    for eachsource in onlyfiles:
        #print(eachsource)
        spectrum = Spectrum(eachsource, 'guo_cmu', calibration_base)
        spectrum.calibrate()
        spectrum.plot(outpath = processedpath,showflag = plotflag)
        spectrum.export(outpath = processedpath)
        
if __name__ == '__main__':
    import sys
    main(sys.argv)
    
    

    
        




