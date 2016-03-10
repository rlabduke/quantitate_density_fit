import os,sys
import argparse
import density_fit_utils
#from iotbx import file_reader

def get_args() :
  # define defaults
  default_add_distance = -0.75

  desc = "Decribe script here"
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdb_file', help='A pdb file')
  parser.add_argument('map_coeff_file', help='Map coefficient file')
  hs = 'specify a residue number to evaluate'
  parser.add_argument('-n','--res_num',dest='res_num',type=int,help=hs)
  hs = 'specify a chain to evaluate'
  parser.add_argument('-c','--chain',dest='chain',help=hs)
  hs = 'a positive or negative value to be added to the default vdw distance'
  hs+= 'when making the surface at which we will calculate map values.'
  hs+= ' Default distance = %s' % default_add_distance
  parser.add_argument('-a','--add_distance',dest='add_distance',
                      type=float,help=hs)
  hs = 'specify an out file name for the kinemage.'
  parser.add_argument('-k','--kinemage_out',dest='kin_out',help=hs)
  args = parser.parse_args()

  assert os.path.exists(args.pdb_file)
  assert os.path.exists(args.map_coeff_file)
  if not args.add_distance : args.add_distance = default_add_distance
  return args

def run(args) :
  args = get_args()

  # Get the density map object
  DM = density_fit_utils.DensityMap(args.pdb_file,args.map_coeff_file)

  # Here is how to get a map value at a specified point.
# xyz = (20.112,37.025,24.574)
# print DM.get_map_value(xyz)

# probe
  pdb_hierarchy = density_fit_utils.get_pdb_hierarchy(args.pdb_file)
  if args.kin_out : fle = open(args.kin_out,'w')
  write_head = True
  for chain in pdb_hierarchy.chains():
    if args.chain and chain.id != args.chain : continue
    for residue_group in chain.residue_groups():
      for conformer in residue_group.conformers():
        for residue in conformer.residues():
          if args.res_num and residue.resseq_as_int() != args.res_num: continue
          # make the probe dot surface for the given residue
          probe_dots = density_fit_utils.ProbeDots(
                         args.pdb_file,
                         chain.id,
                         residue.resseq_as_int(),
                         radius_add = args.add_distance)
          # calculate the map values at each dot
          probe_dots.get_xyz_density_values(DM)
          # write kin
          if args.kin_out :
            probe_dots.write_kin_dots_and_density_values(log=fle)
          # write score csv to stdout
          probe_dots.write_comprehencive_score(write_head=write_head)
          if write_head : write_head = False

  if args.kin_out : fle.close()

if __name__ == '__main__' :
  run(sys.argv[1:])

