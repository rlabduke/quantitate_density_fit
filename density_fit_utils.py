import sys
import iotbx.phil
from cctbx import miller
from iotbx import file_reader
import mmtbx.utils
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
import iotbx.pdb
from libtbx import group_args
from cctbx.array_family import flex
from cStringIO import StringIO
from cctbx import maptbx
from iotbx import mtz
import subprocess
import mmtbx.maps

core_params_str = """\
atom_radius = None
  .type = float
  .help = Atomic radius for map CC calculation. Determined automatically if \
          if None is given
  .expert_level = 2
hydrogen_atom_radius = None
  .type = float
  .help = Atomic radius for map CC calculation for H or D.
  .expert_level = 2
resolution_factor = 1./4
  .type = float
use_hydrogens = None
  .type = bool
"""

master_params_str = """\
%s
scattering_table = *n_gaussian wk1995 it1992 neutron
  .type = choice
  .help = Scattering table for structure factors calculations
detail = atom residue *automatic
  .type = choice(multi=False)
  .help = Level of details to show CC for
map_1
  .help = First map to use in map CC calculation
{
 type = Fc
    .type = str
    .help = Electron density map type. Example xmFobs-yDFcalc (for \
            maximum-likelihood weighted map) or xFobs-yFcalc (for simple \
            unweighted map), x and y are any real numbers.
  fill_missing_reflections = False
    .type = bool
  isotropize = False
    .type = bool
  grid_resolution_factor = 1/4.
}
map_2
  .help = Second map to use in map CC calculation
{
 type = 2mFo-DFc
   .type = str
   .help = Electron density map type. Example xmFobs-yDFcalc (for \
           maximum-likelihood weighted map) or xFobs-yFcalc (for simple \
           unweighted map), x and y are any real numbers.
 fill_missing_reflections = True
   .type = bool
 isotropize = True
   .type = bool
  grid_resolution_factor = 1/4.
}
pdb_file_name = None
  .type = str
  .help = PDB file name.
reflection_file_name = None
  .type = str
  .help = File with experimental data (most of formats: CNS, SHELX, MTZ, etc).
data_labels = None
  .type = str
  .help = Labels for experimental data.
high_resolution = None
  .type=float
low_resolution = None
  .type=float
pdb_id  = None
  .type = str
  .help = the pdb id for csv output
show_csv = False
   .type = bool
   .help = write csv to stdout
show_default_output = False
   .type = bool
"""%core_params_str

def master_params():
  return iotbx.phil.parse(master_params_str, process_includes=False)

def broadcast(m, log = sys.stdout) :
  print >> log, "-"*79
  print >> log, m
  print >> log, "*"*len(m)

def extract_data_and_flags(params, crystal_symmetry=None):
  data_and_flags = None
  if(params.reflection_file_name is not None):
    reflection_file = reflection_file_reader.any_reflection_file(
      file_name = params.reflection_file_name)
    reflection_file_server = reflection_file_utils.reflection_file_server(
      crystal_symmetry = crystal_symmetry,
      force_symmetry   = True,
      reflection_files = [reflection_file])
    parameters = mmtbx.utils.data_and_flags_master_params().extract()
    parameters.force_anomalous_flag_to_be_equal_to = False
    if(params.data_labels is not None):
      parameters.labels = [params.data_labels]
    #### Not Relevant
    # if(params.high_resolution is not None):
      # parameters.high_resolution = params.high_resolution
    # if(params.low_resolution is not None):
      # parameters.low_resolution = params.low_resolution
    #### Not Relevant #### 
    data_and_flags = mmtbx.utils.determine_data_and_flags(
      reflection_file_server = reflection_file_server,
      parameters             = parameters,
      data_description       = "X-ray data",
      extract_r_free_flags   = False, # XXX
      log                    = StringIO())
  return data_and_flags

def pdb_to_xrs(pdb_file_name, scattering_table='n_gaussian'):
  pdb_inp = iotbx.pdb.input(file_name = pdb_file_name)
  xray_structure = pdb_inp.xray_structure_simple()
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq() # VERY important to do.
  mmtbx.utils.setup_scattering_dictionaries(
    scattering_table = scattering_table,
    xray_structure = xray_structure,
    d_min = None)
  return group_args(
    xray_structure = xray_structure,
    pdb_hierarchy  = pdb_hierarchy)

def get_data_labels(mtz_file) :
  mtz_object = mtz.object(file_name=mtz_file)
  labels = mtz_object.column_labels()
  if 'FOBS' in labels and 'SIGFOBS' in labels : return 'FOBS,SIGFOBS'

def get_fmodel_pdb_hierarchy(pdb_file,mtz_file,params,log=sys.stderr) :
  params.pdb_file_name = pdb_file
  params.reflection_file_name = mtz_file
  params.data_labels = get_data_labels(mtz_file)
# broadcast(m='Data labels : %s' % params.data_labels)
  broadcast(m="Input PDB file name: %s"%params.pdb_file_name, log=log)
  pdbo = pdb_to_xrs(pdb_file_name=params.pdb_file_name,
    scattering_table=params.scattering_table)
  pdbo.xray_structure.show_summary(f=log, prefix="  ")
  m = "Input reflection file name: %s"
  broadcast(m=m % params.reflection_file_name, log=log)
  data_and_flags = extract_data_and_flags(params = params)
  data_and_flags.f_obs.show_comprehensive_summary(f=log, prefix="  ")
  xray_structure = pdbo.xray_structure
  # remove hydrogens
  pdb_hierarchy = get_pdb_hierarchy(params.pdb_file_name)
  sele = ~(xray_structure.hd_selection())
  n_hd = sele.count(False)
  if (n_hd > 0) :
    pdb_hierarchy = pdb_hierarchy.select(sele)
    xray_structure = xray_structure.select(sele)
  # create fmodel
  r_free_flags = data_and_flags.f_obs.array(
    data = flex.bool(data_and_flags.f_obs.size(), False))
  fmodel = mmtbx.utils.fmodel_simple(
    xray_structures     = [xray_structure],
    scattering_table    = params.scattering_table,
    f_obs               = data_and_flags.f_obs,
    r_free_flags        = r_free_flags)
  return fmodel,pdb_hierarchy

def set_radius(d_min):
  if(d_min < 1.0):                    atom_radius = 1.0
  elif(d_min >= 1.0 and d_min<2.0):   atom_radius = 1.5
  elif(d_min >= 2.0 and d_min < 4.0): atom_radius = 2.0
  else:                               atom_radius = 2.5
  return atom_radius

def get_pdb_hierarchy(pdb_file) :
  pdb_in = file_reader.any_file(pdb_file)
  pdb_in.check_file_type('pdb')
  input_pdb = iotbx.pdb.input(file_name=pdb_file)
  return input_pdb.construct_hierarchy()

class DensityMap(object) :
  def __init__(self,pdb_file,reflection_file,log=sys.stdout) :
    self.pdb_file = pdb_file
    self.reflection_file = reflection_file
    self.set_map()

  def set_map(self) :
    params = master_params().extract()
    self.fmodel,self.pdb_hierarchy = get_fmodel_pdb_hierarchy(
                                  pdb_file   = self.pdb_file,
                                  mtz_file   = self.reflection_file,
                                  params     = params)
    self.unit_cell = self.fmodel.xray_structure.unit_cell()
    e_map_obj = self.fmodel.electron_density_map()
    coeffs = e_map_obj.map_coefficients(
      map_type     = '2mFo-DFc',
      fill_missing = True,
      isotropize   = True)
    fft_map = coeffs.fft_map(resolution_factor = params.resolution_factor)
    fft_map.apply_sigma_scaling()
    self.map = fft_map.real_map_unpadded()

  def get_map_value(self,xyz) :
    return self.map.eight_point_interpolation(
      self.unit_cell.fractionalize(xyz))

class ProbeDots(object) :
  def __init__(self,file_name,chain,resseq,altid=None,radius_add=-0.75) :
    self.file_name  = file_name
    self.chain      = chain
    self.resseq     = resseq
    self.altid      = altid
    self.radius_add = radius_add
    self.xyzs       = []
    self.run_probe_and_deposit_xyzs()

  def run_probe_and_deposit_xyzs(self) :
    probesele = 'CHAIN_%s %i' % (self.chain,self.resseq)
    args = ['phenix.probe','-q','-drop','-rad0.0','-add%.2f' % self.radius_add]
    #args = ['phenix.probe','-q','-drop','-rad0.0']
    args+= ['-out', probesele, self.file_name]
    pop = subprocess.Popen(args,stdout=subprocess.PIPE)
    self.kinemage_str = pop.communicate()[0]
    assert self.kinemage_str.strip() != ''
    self.deposit_xyzs()

  def deposit_xyzs(self) :
    there = False
    for l in self.kinemage_str.split('\n') :
      if not there or l.strip() == '' :
        if l.startswith('@dotlist') : there = True
        continue
      if not l.startswith('{') : continue
      sl = l[l.rfind(' '):].strip()
      xyzstr = sl.split(',')
      assert len(xyzstr) == 3, xyzstr
      xyz = (float(xyzstr[0]),float(xyzstr[1]),float(xyzstr[2]))
      self.xyzs.append(xyz)

  def get_xyz_density_values(self,DensityMap) :
    self.xyz_density_values = []
    for xyz in self.xyzs :
      v = DensityMap.get_map_value(xyz)
      self.xyz_density_values.append((v,xyz))
      #print v,xyz
    self.set_density_fit_score()

  def set_density_fit_score(self) :
    self.density_fit_score = 0
    self.dots_n = 0
    self.dots_lt1, self.dots_gte1 = 0,0
    for v,xyz in self.xyz_density_values :
      self.dots_n += 1
      if  v < 1.0 :
        self.dots_lt1 += 1
        self.density_fit_score += v
      else : self.dots_gte1 += 1

  def write_kin_dots(self,log=sys.stdout) :
    print >> log, self.kinemage_str

  def write_kin_dots_and_density_values(self,log=sys.stdout) :
    col = 'gray'
    print >> log, '''@group dominant {%s,%s,%s}
@subgroup dominant {extern dots}
@master {surface}
@balllist {x} color=white radius= 0.02 master={surface} nohilite''' % (self.chain,
                                                self.resseq,
                                                self.altid)
    line = '{2fo-fc = %.4f} %s %.3f,%.3f,%.3f'
    for v,xyz in self.xyz_density_values :
      if   v < 0.8 : col = 'red'
      elif v < 0.9 : col = 'orange'
      elif v < 1.0 : col = 'yellow'
      elif v < 1.1 : col = 'green'
      else         : col = 'blue'
      print >> log, line % (v,col,xyz[0],xyz[1],xyz[2])

  def write_comprehencive_score(self,write_head=False,log=sys.stdout) :
    head = "chain,resseq,altid,dots_n,dots_lt1,dots_gte1,score"
    if write_head : print >> log, head
    print >> log, '%s,%i,%s,%i,%i,%i,%.4f' % (self.chain,self.resseq,self.altid,
                           self.dots_n, self.dots_lt1, 
                           self.dots_gte1, self.density_fit_score)
