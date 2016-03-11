import os,sys
import argparse
import density_fit_utils
#from iotbx import file_reader
import pickle

def get_args() :
  # define defaults
  default_vdw_scale = 0.5
  default_map_scale = 'sigma'
  default_map_type = '2mFo-DFc'

  desc = "Decribe script here"
  parser = argparse.ArgumentParser(description=desc)
  parser.add_argument('pdb_file', help='A pdb file')
  parser.add_argument('map_coeff_file', help='Map coefficient file')
  hs = 'specify a residue number to evaluate'
  parser.add_argument('-n','--res_num',dest='res_num',type=int,help=hs)
  hs = 'specify a chain to evaluate'
  parser.add_argument('-c','--chain',dest='chain',help=hs)
  hs = 'a scale factor for the default vdw distance to determine the distance'
  hs+= 'of the new surface where map values will bel calculated.'
  hs+= ' Default scale = %s' % default_vdw_scale
  parser.add_argument('-a','--vdw_scale',dest='vdw_scale',
                      type=float,help=hs)
  hs = 'specify an out file name for the kinemage.'
  parser.add_argument('-k','--kinemage_out',dest='kin_out',help=hs)
  hs = 'create shells of surfaces around each residue'
  parser.add_argument('-s','--shells',dest='shells',help=hs,action='store_true')
  hs = 'specify map scale, either sigma or volume - sigma is deault'
  parser.add_argument('-x','--map_scale',dest='map_scale',help=hs)
  hs = 'specify a map type'
  parser.add_argument('-m','--map_type',dest='map_type',help=hs)
  args = parser.parse_args()

  assert os.path.exists(args.pdb_file)
  assert os.path.exists(args.map_coeff_file)
  if not args.vdw_scale : args.vdw_scale = default_vdw_scale
  if not args.map_scale : args.map_scale = default_map_scale
  if not args.map_type : args.map_type   = default_map_type
  return args

def run(args) :
  args = get_args()

  # Get the density map object
  base = os.path.basename(args.pdb_file)
  picklefn = base[:base.rfind('.')] + '.pickle'
  if os.path.exists(picklefn) :
    fle = open(picklefn,'r')
    DM = pickle.load(fle)
    fle.close()
  else :
    DM = density_fit_utils.DensityMap(pdb_file = args.pdb_file,
                                      reflection_file = args.map_coeff_file,
                                      map_scale = args.map_scale,
                                      map_type  = args.map_type)
    fle = open(picklefn,'w')
    pickle.dump(DM,fle)
    fle.close()

  # Here is how to get a map value at a specified point.
# xyz = (20.112,37.025,24.574)
# print DM.get_map_value(xyz)

# probe
  pdbo = density_fit_utils.pdb_to_xrs(args.pdb_file)
  pdb_hierarchy  = pdbo.pdb_hierarchy
  xray_structure = pdbo.xray_structure
  # eliminate hydrogens
  sele = ~(xray_structure.hd_selection())
  n_hd = sele.count(False)
  if (n_hd > 0) :
    pdb_hierarchy = pdb_hierarchy.select(sele)
    xray_structure = xray_structure.select(sele)
  nohfn = base[:base.rfind('.')] + '_noH' + base[base.rfind('.'):]
  pdb_hierarchy.write_pdb_file(nohfn)
  scalelist = [args.vdw_scale]
  if args.shells : scalelist = [0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
  if args.kin_out : fle = open(args.kin_out,'w')
  write_head = True
  for chain in pdb_hierarchy.chains():
    if args.chain and chain.id != args.chain : continue
    for residue_group in chain.residue_groups():
      for conformer in residue_group.conformers():
        for residue in conformer.residues():
          if args.res_num and residue.resseq_as_int() != args.res_num: continue
          # make the probe dot surface for the given residue
          for scale in scalelist :
            probe_dots = density_fit_utils.ProbeDots(
                           nohfn,
                           chain.id,
                           residue.resseq_as_int(),
                           radius_scale = scale)
            # calculate the map values at each dot
            probe_dots.get_xyz_density_values(DM)
            # write kin
            if args.kin_out :
              probe_dots.write_kin_dots_and_density_values(log=fle)
            # write score csv to stdout
            probe_dots.write_comprehencive_score(write_head=write_head,
                                                 format='human')
            if write_head : write_head = False

  if args.kin_out : fle.close()

if __name__ == '__main__' :
  run(sys.argv[1:])

