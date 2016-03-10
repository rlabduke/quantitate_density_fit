# quantitate_density_fit
This is under heave developement as of 03-10-2016

Note that you must source PHENIX in order to run this:
  $ source /path/to/phenix/build/setpaths.sh

To see help on running run :
  $ python quantitavie_ED_fit.py -h

You need a pdb and structure factor mtz to run. These can be fetched and calculated via phenix :
  $ phenix.fetch_pdb --mtz xxxx

where xxxx is a pdb code.

Notes that if you are making a kinemage, it is a good idea to specify a chain (-c) and resedue sequence number (-n) else you get a VERY large kinemage file.
