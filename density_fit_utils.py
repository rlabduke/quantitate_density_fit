import sys,os
from iotbx.file_reader import any_file
import mmtbx.utils
import iotbx.pdb
from cctbx.array_family import flex
from iotbx import mtz
import matplotlib.pyplot as plt
from cStringIO import StringIO
import subprocess
import pickle

def get_data_labels(mtz_file) :
  mtz_object = mtz.object(file_name=mtz_file)
  labels = mtz_object.column_labels()
  if 'FOBS' in labels and 'SIGFOBS' in labels : return 'FOBS,SIGFOBS'

def write_line_plot(x,y,xlabel,ylabel,filename,ylim=None,
                    title=None,subtitle=None) :
  fig = plt.figure() # needed if you want to save the figure
  # Fake labels
  if subtitle : fig.suptitle(subtitle,fontsize=18)
  plt.plot(x,y)
  plt.xlabel(xlabel)
  plt.ylabel(ylabel)
  if ylim :
    assert len(ylim) == 2
    assert type(ylim) == tuple
    plt.ylim(ylim)
  if not title : title = '%s vs. %s' % (xlabel,ylabel)
  plt.title(title)
  if not filename.endswith('.png') : filename = filename + '.png'
  fig.savefig(filename, format='png')
# plt.show()# this bring up a gui
  plt.close()

class DensityMap(object) :
  def __init__(self,pdb_file,reflection_file,
               map_scale='sigma',map_type='2mFo-DFc',log=sys.stdout) :
    assert map_scale in ['volume','sigma']
    self.map_scale = map_scale
    self.map_type = map_type
    self.pdb_file = pdb_file
    self.log = log
    self.reflection_file = reflection_file
    self.set_map()

  def set_map(self) :
    pdb_in = any_file(self.pdb_file)#W
    self.pdb_hierarchy = pdb_in.file_object.hierarchy#W
    inputs = mmtbx.utils.process_command_line_args([self.pdb_file,
                                                  self.reflection_file])
    data_labels = get_data_labels(self.reflection_file)
    parameters = mmtbx.utils.data_and_flags_master_params().extract()
    parameters.force_anomalous_flag_to_be_equal_to = False
    parameters.labels = [data_labels]

    data_and_flags = mmtbx.utils.determine_data_and_flags(
      reflection_file_server = inputs.get_reflection_file_server(),
      parameters             = parameters,
      extract_r_free_flags   = False, # XXX
      keep_going             = True) # don't stop if free flags are not present
    f_obs = data_and_flags.f_obs
    r_free_flags = data_and_flags.f_obs.array(
          data = flex.bool(data_and_flags.f_obs.size(), False))
    xrs = iotbx.pdb.input(
      file_name=inputs.pdb_file_names[0]).xray_structure_simple()
    self.fmodel = mmtbx.utils.fmodel_simple(
      f_obs=f_obs,
      r_free_flags=r_free_flags,
      scattering_table="n_gaussian",
      xray_structures=[xrs])
    self.unit_cell = self.fmodel.xray_structure.unit_cell() 
    ###
    self.unit_cell = self.fmodel.xray_structure.unit_cell()
    e_map_obj = self.fmodel.electron_density_map()
    coeffs = e_map_obj.map_coefficients(
      map_type     = '2mFo-DFc',
      fill_missing = True,
      isotropize   = True)
    fft_map = coeffs.fft_map(resolution_factor = 1./4)
    fft_map.apply_sigma_scaling()
    self.map = fft_map.real_map_unpadded()

  def get_map_value(self,xyz) :
    return self.map.eight_point_interpolation(
      self.unit_cell.fractionalize(xyz))

class ProbeLabel(object) :
  def __init__(self,label_str) :
    self.label_str = label_str
    self.atom = label_str[:4]
    self.alt_loc = label_str[4]
    self.res_type = label_str[5:8]
    self.res_num = int(label_str[8:12].strip())
    self.icode = label_str[12]
    self.chain = label_str[13:]

  def __str__(self) :
    return self.label_str

class ProbeDots(object) :
  def __init__(self,
               file_name,
               chain,
               resseq,
               altloc=None,
               icode=None,
               radius_scale=-0.75,
               include_waters=False) :
    self.file_name      = file_name
    self.chain          = chain
    self.resseq         = resseq
    self.altloc         = altloc
    self.icode          = icode
    self.radius_scale   = radius_scale
    self.include_waters = include_waters
    self.xyzs           = []
    self.xyz_density_values = {}
    self.bb_atoms = [' N  ',' CA ',' O  ',' C  ']
    # analytics dict, contains all
    self.density_fit_score = {}
    self.dots_n = {}
    self.dots_lt1 = {}
    self.dots_gte1 = {}
    self.average_density = {}
    self.dots_lt1_ratio = {}
    self.dots_gte1_ratio = {}
    self.run_probe_and_deposit_xyzs()

  def run_probe_and_deposit_xyzs(self) :
    probesele = 'CHAIN_%s %i' % (self.chain,self.resseq)
    probesele = 'CHAIN_%s %i-%i' % (self.chain,self.resseq-1,self.resseq+1)
    args = ['phenix.probe','-q','-rad0.0','-scale%.2f' % self.radius_scale]
    if not self.include_waters : args += ['-nowat']
    #args = ['phenix.probe','-q','-drop','-rad0.0']
    args+= ['-drop','-out', probesele, self.file_name]
    print >> sys.stderr, ' '.join(args)
    pop = subprocess.Popen(args,stdout=subprocess.PIPE)
    self.kinemage_str = pop.communicate()[0]
    assert self.kinemage_str.strip() != ''
    self.deposit_xyzs()

  def deposit_xyzs(self) :
    there = False
    Label = None
    for l in self.kinemage_str.split('\n') :
      if not there or l.strip() == '' :
        if l.startswith('@dotlist') : there = True
        continue
      if not l.startswith('{') : continue
      if not l.startswith('{"}') : label = ProbeLabel(l[1:l.find('}')])
      assert label != None
      if label.res_num != self.resseq : continue
      sl = l[l.rfind(' '):].strip()
      xyzstr = sl.split(',')
      assert len(xyzstr) == 3, xyzstr
      xyz = ((float(xyzstr[0]),float(xyzstr[1]),float(xyzstr[2])),label)
      self.xyzs.append(xyz)

  def get_xyz_density_values(self,DensityMap) :
    mt = DensityMap.map_type
    self.xyz_density_values[mt] = {}
    for xyz,label in self.xyzs :
      v = DensityMap.get_map_value(xyz)
      self.xyz_density_values[mt][xyz] = (v,label)
      #print v,xyz
    self.set_density_fit_score(typ='all')
    self.set_density_fit_score(typ='sc')
    self.set_density_fit_score(typ='bb')

  def set_density_fit_score(self,typ) :
    # typ can be bb, sc, or all
    addv = 0
    self.density_fit_score[typ] = 0
    self.dots_n[typ] = 0
    self.dots_lt1[typ], self.dots_gte1[typ] = 0,0
    for xyz,vl in self.xyz_density_values['2mFo-DFc'].items() :
      v,label = vl
      if typ == 'sc' and label.atom in self.bb_atoms : continue
      if typ == 'bb' and label.atom not in self.bb_atoms : continue
      addv += v
      self.dots_n[typ] += 1
      if  v < 1.0 :
        self.dots_lt1[typ] += 1
        self.density_fit_score[typ] += v
      else : self.dots_gte1[typ] += 1
    self.average_density[typ] = addv/self.dots_n[typ]
    self.dots_lt1_ratio[typ] = (self.dots_lt1[typ]*1.0)/self.dots_n[typ]
    self.dots_gte1_ratio[typ] = (self.dots_gte1[typ]*1.0)/self.dots_n[typ]

  def write_kin_dots(self,log=sys.stdout) :
    print >> log, self.kinemage_str

  def write_kin_dots_and_density_values(self,log=sys.stdout) :
    col = 'gray'
    print >> log, '''@group {%s,%s,%s,%.2f} animate dominant
@subgroup dominant {extern dots}
@master {surface}
@pointmaster 'S' {sidechain}
@pointmaster 'B' {backbone}
@balllist {x} color=white radius= 0.02 master={surface} nohilite''' % (
                                                self.chain,
                                                self.resseq,
                                                self.altloc,
                                                self.radius_scale)
    line = '{%s 2fo-fc = %.4f} \'%s\' %s %.3f,%.3f,%.3f'
    for xyz,vl in self.xyz_density_values['2mFo-DFc'].items() :
      v,label = vl
      if label.atom in self.bb_atoms : scbb = 'B'
      else : scbb = 'S'
      if   v < 0.5 : col = 'red'
      elif v < 0.6 : col = 'orange'
      elif v < 0.7 : col = 'gold'
      elif v < 0.8 : col = 'yellow'
      elif v < 0.9 : col = 'lime'
      elif v < 1.0 : col = 'green'
      elif v < 1.1 : col = 'cyan'
      elif v < 1.2 : col = 'sky'
      elif v < 1.3 : col = 'blue'
      else         : col = 'purple'
      print >> log, line % (label.__str__(),v,scbb,col,xyz[0],xyz[1],xyz[2])

  def write_comprehensive_score(self,write_head=False,
                                format='csv',typ='all',log=sys.stdout) :
    assert format in ['csv','human']
    head = "chain,resseq,altloc,surf_scale,dots_n,dots_lt1,dots_gte1,score,"
    head+= "avg_density"
    llist = [self.chain,self.resseq,self.altloc,self.radius_scale]
    llist+= [self.dots_n['all'],self.dots_lt1['all'],self.dots_gte1['all']]
    llist+= [self.density_fit_score['all'],self.average_density['all']]
    csvs = '%s,%i,%s,%.2f,%i,%i,%i,%.4f,%.4f' % tuple(llist)
    if format == 'csv' :
      if write_head : print >> log, head
      print >> log, csvs
    else :
      heads = head.split(',')
      l = csvs.split(',')
      ls = []#e.ljust(10) for e in l]
      headformat = []
      assert len(l) == len(heads)
      for i,e in enumerate(l) :
        if heads[i] == 'score' : ljl = 14
        elif len(heads[i]) < 10 : ljl = 10
        else : ljl =len(heads[i])
        ls.append(e.ljust(ljl))
        headformat.append(heads[i].ljust(ljl))
      if write_head : print >> log, ''.join(headformat)
      print >> log, ''.join(ls)

class ResidueDensityShells(object) :

  def __init__(self,pdb_file,mtz_file,map_scale,
               chain=None,res_num=None,alt_loc=None,
               analysis_type='2mFo-DFc') :
    self.pdb_file = pdb_file
    self.basename = os.path.basename(self.pdb_file)
    self.mtz_file = mtz_file
    self.map_scale = map_scale
    self.chain = chain
    self.res_num = res_num
    self.alt_loc = alt_loc
    self.analysis_type = analysis_type
    self.scalelist = [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1]
    assert self.analysis_type in ['2mFo-DFc','correlation','comprehensive']
    self.DM_tfo_fc = self.get_map(map_type='2mFo-DFc')
    if analysis_type != '2mFo-DFc' : self.DM_fc = self.get_map(map_type='Fc')
    self.set_hierarchy()
    self.set_density_shells()

  def get_map(self,map_type='2mFo-DFc') :
    print >> sys.stderr, 'calculating %s map' % map_type
    picklefn = self.basename[:self.basename.rfind('.')] + '_%s.pickle'%map_type
    if os.path.exists(picklefn) :
      fle = open(picklefn,'r')
      DM = pickle.load(fle)
      fle.close()
    else :
      DM = DensityMap(pdb_file = self.pdb_file,
                                        reflection_file = self.mtz_file,
                                        map_scale = self.map_scale,
                                        map_type  = map_type,
                                        log       = StringIO())
      fle = open(picklefn,'w')
      DM.log = None# Needs to be done because you cannot pickle file objects
      pickle.dump(DM,fle)
      fle.close()
    return DM

  def set_hierarchy(self) :
    pdb_in = any_file(self.pdb_file)
    self.pdb_hierarchy = pdb_in.file_object.hierarchy
    xray_structure = iotbx.pdb.input(
      file_name=self.pdb_file).xray_structure_simple()
    # eliminate hydrogens
    sele = ~(xray_structure.hd_selection())
    n_hd = sele.count(False)
    if (n_hd > 0) :
      self.pdb_hierarchy = self.pdb_hierarchy.select(sele)
      xray_structure = xray_structure.select(sele)
    self.nohfn = self.basename[:self.basename.rfind('.')] + '_noH'
    self.nohfn += self.basename[self.basename.rfind('.'):]
    self.pdb_hierarchy.write_pdb_file(self.nohfn)

  def set_density_shells(self) :
    self.density_shells = []
    for chain in self.pdb_hierarchy.chains():
      if self.chain and chain.id != self.chain : continue
      for residue_group in chain.residue_groups():
        for conformer in residue_group.conformers():
          for residue in conformer.residues():
            if self.res_num and residue.resseq_as_int() != self.res_num:continue
            if residue.resname == 'HOH' : include_waters = True
            else : include_waters = False
            # make the probe dot surface for the given residue
            for scale in self.scalelist :
              probe_dots = ProbeDots(
                           self.nohfn,
                           chain.id,
                           residue.resseq_as_int(),
                           altloc=conformer.altloc,
                           icode=residue.icode,
                           radius_scale = scale,
                           include_waters = include_waters)
	      probe_dots.get_xyz_density_values(self.DM_tfo_fc)
              if self.analysis_type != '2mFo-DFc' :
                probe_dots.get_xyz_density_values(self.DM_fc)
              self.density_shells.append((scale,probe_dots))

  def write_shells_kin(self,log=sys.stdout) :
    for scale,probe_dots in self.density_shells :
      probe_dots.write_kin_dots_and_density_values(log=log)

  def write_shells_scores(self,log=sys.stdout) :
    write_head = True
    for scale,probe_dots in self.density_shells :
      probe_dots.write_comprehensive_score(write_head=write_head,
                                                 format='human')
      if write_head : write_head = False

  def write_plots(self,prefix,typ) :
    # typ can be bb, sc, or all
    val_dict = {'radius_scale':[],'dots_n':[],'dots_lt1':[],'dots_gte1':[]}
    val_dict['dots_lt1_ratio'] = []
    val_dict['dots_gte1_ratio'] = []
    val_dict['density_fit_score'] = []
    val_dict['average_density'] = []
    for scale,probe_dots in self.density_shells :
      val_dict['radius_scale'].append(probe_dots.radius_scale)
      val_dict['dots_n'].append(probe_dots.dots_n[typ])
      val_dict['dots_lt1'].append(probe_dots.dots_lt1[typ])
      val_dict['dots_lt1_ratio'].append(probe_dots.dots_lt1_ratio[typ])
      val_dict['dots_gte1'].append(probe_dots.dots_gte1[typ])
      val_dict['dots_gte1_ratio'].append(probe_dots.dots_gte1_ratio[typ])
      val_dict['density_fit_score'].append(probe_dots.density_fit_score[typ])
      val_dict['average_density'].append(probe_dots.average_density[typ])
    tn = ['dots_n','dots_lt1','dots_lt1_ratio','dots_gte1','dots_gte1_ratio']
    tn +=['density_fit_score','average_density']
    for type_number in tn :
      if type_number == 'average_density' : ylim = (-1,6)
      else : ylim = None
      fn = "%s_%s_%s" % (prefix,type_number,typ)
      title = '%s vs. %s %s' % ('radius_scale',type_number,typ)
      write_line_plot(x=val_dict['radius_scale'],
                      y=val_dict[type_number],
                      xlabel='radius_scale',
                      ylabel=type_number,
                      title = title,
                      filename=fn,
                      ylim=ylim)
