import pandas as pd

def GetMarkerList(configFile, selectedColumn):
################################################
# From tsv file containing markers and flags
# return a list of marker's that are flagged as 1
# input
# file, column name
#
# output:
# list of matched markers
################################################
    markerColumn = "marker_name"
    df = pd.read_csv(configFile, sep="\t")
    if (selectedColumn in df):
        filtered = (df[df[selectedColumn] == 1.0])
        markers = filtered[markerColumn].values.tolist()
        return(markers)
    else:
        print(selectedColumn + " does not exist as a column heading!")
        return[] 

#a = GetMarkerList("/t1-data/project/covidhyperion/staylor/covid2/config/markers_full.tsv","nucleus")
#print(a)