# mBUD

Metal-organic framework (MOF) has relevance in extensive applications such as gas adsorption, separation, and energy storage. It constitutes metal nodes, organic linkers, and functional groups as its building units (BUs). 
Here, we present the MOF Building Unit Developer (mBUD) package that extracts the BUs from a given MOF structure. The code for this package is available on our site. Along with this, the web-version of mBUD has also been provided.
We have employed mBUD to extract the BUs from the Computation-Ready, Experimental (CoRE) MOF 2019 database. The BUs database obtained from this analysis is provided for further utilization, such as application-specific MOF development.

## Code

We have provided the code for the entire package here as well as on our [web page](https://cnislab.com/mbud/code). This is a Python version 3.7 based object-oriented code. It consists of three prime classes library, atom, and MOF. The library incorporates the input data required for the analysis (such as atom typing radius (ATR), skin distance, etc.). The atom class includes all the relevant information and functions for the atoms, such as the list of neighbor atoms, metal/non-metal classification, etc. Likewise, the MOF class contains the functions applicable to the MOF structure like solvent identification, bonded network formation, deconstruction of MOF, etc. The final output is the extracted building units in the Crystallographic Information File (CIF) format.

__First install the dependencies via conda or pip:__
* pandas>=1.2.4
* numpy>=1.20.1
* networkx>=2.5
* datetime

__Use the following command to run the code:__
```
python3 mBUD.py --ciffile path/ABCDEF.cif --outdir path/output
```
All the extracted building units will be stored in the _output_ folder.

## Database

The building unit (BU) [database](https://cnislab.com/mbud/database) consists of all the unique BU extracted from the subset Computation-Ready, Experimental (CoRE) MOF 2019-ASR database. The subset of CoRE database consists of 9,268 MOF structures. We extracted  2,580 BUs ( including metal nodes, and organic linkers) from this subset. 
The provided database comprises both the experimental and computational BUs. Experimental BUs are essential for the visualization of the MOF chemistry. On the other hand, the computational BUs can be readily employed to construct MOF crystals computationally.

## Web

We have also presented the web-version of the [mBUD package](https://cnislab.com/mbud/tool). To extract the building units (BUs) from a MOF structure, upload the Crystallographic Information File (CIF) for the MOF. Once the analysis is complete, the extracted BUs (in both experimental and computational representation) can be visualized and downloaded in the CIF format.
Currently, our code can read 'P 1' space group only. So, the Crystallographic Information File (CIF) for all the structures need to be converted into the 'P 1' space group. This operation can be accomplished using the cif2cell python package.


