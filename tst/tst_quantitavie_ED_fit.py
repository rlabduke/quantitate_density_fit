import os,sys
sys.path.append('..')
import density_fit_utils
import cStringIO

def run(args) :
  pdb_file = '4rp6.pdb'
  hkl_file = '4rp6.mtz'
  chain = 'Z'
  res_num = 2
  rds = density_fit_utils.ResidueDensityShells(pdb_file = pdb_file,
                                               mtz_file = hkl_file,
                                               map_scale = 'sigma',
                                               chain = chain,
                                               res_num = res_num,
                                               analysis_type = 'correlation')

  files = [f for f in os.listdir('.') if os.path.isfile(f)]
  for f in files:
    if f.endswith('.pickle') :
      os.remove(f)
  prefix = "%s_%s_%s" % (pdb_file[:4],chain,res_num)
  fn = "%s_shells.kin" % prefix
  sio = cStringIO.StringIO()
  rds.write_shells_kin(log=sio)
  if not os.path.exists(fn) :
    fle = open(fn,'w')
    print >> fle, sio.getvalue()
    fle.close()
    print '%s written.' % fn
  else :
    failmsg = "Your test failed! "
    fle = open(fn,'r')
    reference_lines = fle.readlines()
    fle.close()
    generated_lines = sio.getvalue().split('\n')
    fm = failmsg+"The kins had different number of lines."
    assert len(reference_lines) == len(generated_lines),fm
    for i in range(len(generated_lines)) :
      if reference_lines[i].strip() != generated_lines[i].strip() :
        fm =  failmsg + "The following lines were different:\n"
        fm +=  'Ref : %s' % reference_lines[i].strip() + '\n'
        fm +=  'Gen : %s' % generated_lines[i].strip() + '\n'
        raise RuntimeError(fm)
    print "\nYour test scceeded! Go have a scotch, preferably one from Islay.\n"

if __name__ == '__main__' :
  run(sys.argv[1:])

