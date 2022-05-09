## mlvupload.py

The script takes a single tab delimited file pf all cells and markers, as well as other metadata and uploads it MLV.  The only argument required is the path to the yaml file which contains all the arguments required. A token, obtained from the MLV website is required for authentocation (see below)


### Basics

To upload data from the pipeline, a lot of parameters in the default yaml file need not br chabge. However, the following parameters require the appropriate values.


* _name_ - the name of the project
* _description_ - brief description of the project
* _genome_ the genome build e.g. hg19, mm10 - not actually used but required
* _token_ - token obtained from the server. In order to  do this, go to the the mlv instance and login. Then navigate to 'My Projects' (on the right of the top Nav Bar).  Cick on the Hyperion project and presse 'get token'. Copy tis token and use it the yaml file. 
* _data_file_ - the path to clusters.txt file generated from the Rpehenoclustering.R script
* _marker_file_ - path to the marker file which describes the markers in the experiment and those used for the initial clustering
* _pca_file_ - the path to the pca.txt file generated from the Rpehenoclustering.R script
* _images_folders_ - a  list of length 1 containing the path to the zegami images produced by the pipeline _.../zegami/img_ 

Then run

    python mlvupload.py /path/to/myconfig.yaml

The data will be compressed and sent to the server where a project will be created. You will receive an email with a link when the project has been built  


