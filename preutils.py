import re
from datetime import datetime
import argparse

#data section

dict_05_22_99 = {
'3e': [0.034683, 121.719, 254.85],
'4e': [0.046211, 123.347, 255.25],
'5e': [0.057866, 124.467, 254.75],
'6e': [0.069429, 125.105, 254.35],
'7e': [0.080945, 125.614, 254.25],
'8e': [0.092428, 125.950, 254.15],
'9e': [0.103754, 126.195, 254.10],
'2e': [0.023750, 113.000, 255.10],
'3r': [0.031136, 124.879, 257.15],
'4r': [0.041499, 125.796, 257.25],
'5r': [0.051936, 126.362, 257.15],
'6r': [0.062293, 126.721, 257.10],
'7r': [0.072651, 126.981, 257.10],
'8r': [0.083060, 127.182, 257.10],
'9r': [0.093264, 127.322, 257.10],
}

dict_06_29_18={
'3e': [0.03231, 121.047, 251.7],
'4e': [0.04315, 122.924, 251.7],
'5e': [0.05394, 124.085, 251.7],
'6e': [0.06506, 124.716, 252.9],
'7e': [0.07572, 125.288, 252.3],
'8e': [0.08635, 125.700, 252.2],
'9e': [0.09713, 126.040, 252.1],
'3r': [0.03109, 124.936, 256.4],
'4r': [0.04152, 125.806, 256.4],
'5r': [0.05178, 126.363, 256.4],
'6r': [0.06207, 126.690, 256.7],
'7r': [0.07238, 126.950, 256.7],
'8r': [0.08275, 127.157, 256.7],
'9r': [0.09303, 127.310, 256.6],
}


calibration_base = {
datetime.strptime('05-22-99', '%m-%d-%y'):dict_05_22_99,
datetime.strptime('06-29-18', '%m-%d-%y'):dict_06_29_18,
}

def munck_parser(source_file):
    rawdataseq = []
    print(source_file)
    with open(source_file,'r') as file_ptr:
        for i, eachline in enumerate(file_ptr):
            if i == 0:
                header1 = eachline.split(' ')
                filename = header1[0]
            if i == 1:
                header2 = eachline
                re_results = re.search('(.+?)/(.+?)/(.+?)\(//\)/(.+?)/(.*)',eachline)
                if re_results is None:
                    re_results = re.search('(.+?)/(.+?)/(.+?)/(.+?)/(.*)',eachline)

                liner = re_results.group(1)
                temperature = re_results.group(2)
                field = re_results.group(3)
                veloscale = re_results.group(4)
                collect_date = re_results.group(5)

            if i == 2:
                pass
            if i == 3:
                pass
                #header4 = eachline.split(' ')
                #print(header4)
                #foldpoint = float(header4[0])
            if i >= 4:
                eight_numbers = eachline.split('.')
                rawdataseq += map(int, eight_numbers[:-1])
        
        collect_date = datetime.strptime(collect_date, '%m-%d-%y')
        wholeblock = filename + ' ' + header2
        
    return filename, collect_date, liner, temperature, field, veloscale, rawdataseq, wholeblock


def read_calibration(calibration_file, collect_date, veloscale):
    #compute the closest date
    daterr = min(calibration_base.keys(), key=lambda dater : abs((collect_date.date()-dater.date())))
    print(veloscale)
    return calibration_base[daterr][veloscale]

def parse_arguments():
    # Command-line flags are defined here.
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-path', dest='sourcepath',
                        type=str, default='./sourcefile',
                        help="Path to the source files' root")
    parser.add_argument('--processed-path', dest='processedpath',
                        type=str, default='./processedfile',
                        help="Path to the processed files' root")
    # ~ parser.add_argument('--calibration', dest='calibration_config_path',
                        # ~ type=str, default='calibration_config_file',
                        # ~ help="Path to the calibration config file.")
    
    parser_group = parser.add_mutually_exclusive_group(required=False)
    parser_group.add_argument('--plot', dest='plotflag',
                              action='store_true',
                              help="Whether to show the plotting promptly")
    parser_group.add_argument('--no-plot', dest='plotflag',
                              action='store_false',
                              help="Whether to show the plotting promptly.")
    parser.set_defaults(render=False)
    return parser.parse_args()

def folddata(foldpoint, channels, dataseq):
    fold_channels = [each if each <=foldpoint else (2*foldpoint - each) for each in channels]
    sorted_pairs = sorted(zip(fold_channels, dataseq), key=lambda x:x[0])
    aved_channel = []
    aved_count = []
    for i in range(len(sorted_pairs)//2):
        this_channel = (sorted_pairs[2*i][0] + sorted_pairs[2*i+1][0])/2
        this_count = (sorted_pairs[2*i][1] + sorted_pairs[2*i+1][1])
        aved_channel.append(this_channel)
        aved_count.append(this_count)
    # ~ from matplotlib import pyplot as plt
    # ~ plt.plot(aved_v, aved_count)
    # ~ plt.show()
    return aved_channel, aved_count
        


            





