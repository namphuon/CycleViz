import sys, random, math
import os
from collections import Set
import bisect
import argparse
import numpy as np
import hg19util as hg19
from collections import defaultdict
import matplotlib
matplotlib.use('Agg')
from decimal import Decimal
import matplotlib.path as mpath
from matplotlib import pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
from matplotlib.collections import PatchCollection
from intervaltree import Interval, IntervalTree


#Create colormap for human chromosomes
chromosome_colors = dict([(hg19.chrName[i],c) for (i,c) in enumerate(plt.cm.get_cmap(None, 23).colors)])
ensembl_grc37, ensembl_grc37_map = load_ensembl_grc37()

contig_spacing = 0.01
seg_spacing = 0.015
bar_width = 2.0/2
global_rot = 90.0
tick_positions = 100000
hg19.load_exons()

outer_bar = 19
bar_width = 2.0/2
bar_spacing = 1
#bed_track_height = 7
bed_track_height = 4
#bed_track_spacing = 1.25
bed_track_spacing = .5
control_delta = -1
minor_tick_delta = 1
gene_delta = 3
font_size = 8

oncogene_map = {}
load_oncogenes()
def load_oncogenes():
  input = open('/pedigree2/projects/namphuon/data/references/hg19/annotations/allOnco_May2018.tsv','r')
  foo = input.next()
  for line in input:
    res = line.replace('"','').split()[1:-1]
    for re in res:
      for r in re.split(','):
        oncogene_map[r.upper()] = res
  input.close()

def load_ensembl_grc37():
  input = open('/pedigree2/projects/namphuon/data/references/hg19/annotations/Homo_sapiens.GRCh37.64.gtf' ,'r')
  ensemble_data = hg19.interval_list()
  ensembl_grc37_map = {}
  for line in input:
    res = line.split('\t')
    info = dict([r.strip().replace('"','').split(' ') for r in res[-1].strip().split('; ') if len(r.split(' ')) == 2])
    temp = hg19.interval("chr%s" % res[0],int(res[3]),int(res[4]),info={'data':info})
    ensemble_data.append(temp)
    foo = ensembl_grc37_map.setdefault(info['gene_id'],[]).append(temp)
  ensembl_grc37 = hg19.interval_list()    
  for g in ensembl_grc37_map:
    start = min([e.start for e in ensembl_grc37_map[g]])
    end = max([e.end for e in ensembl_grc37_map[g]])
    name = [e.info['data']['gene_name'] for e in ensembl_grc37_map[g] if 'gene_name' in e.info['data']]
    name = name[0] if len(name) >= 1 else g          
    ensembl_grc37.append(hg19.interval("%s" % e.chrom, start, end, info={'intervals':ensembl_grc37_map[g],'GeneID':g,'Name':name}))
  input.close()
  ensembl_grc37.sort()
  for e in ensembl_grc37:
    ensembl_grc37_map[e.info['GeneID']] = e
  return (ensembl_grc37, ensembl_grc37_map)


def load_fpkm_file(fpkm_file):
  lines = [line.strip().split('\t') for line in open(fpkm_file, 'r')]
  return dict([(res[1].split('.')[0],res[2]) for res in lines])

def build_genebed_from_fpkm(fpkm):
  fpkm_bed = hg19.interval_list()
  for (g,f) in fpkm.items():
    if g not in ensembl_grc37_map:
      continue
    gene = ensembl_grc37_map[g]
    fpkm_bed.append(hg19.interval(gene.chrom, gene.start, gene.end, info={'value':f, 'name':gene.info['Name']}))
  fpkm_bed.sort()
  return fpkm_bed

def cart2pol(x, y):
  rho = np.sqrt(x**2 + y**2)
  phi = np.arctan2(y, x) / (2. * np.pi) * 360 
  return(rho, phi)

def pol2cart(rho, phi):
  x = rho * np.cos(phi) 
  y = rho * np.sin(phi)
  return(x, y)

def polar_series_to_cartesians(line_points,r):
  x_list = []
  y_list = []
  for i in line_points:
    x,y = pol2cart(r,i)
    x_list.append(x)
    y_list.append(y)
  return x_list,y_list

def rel_genes(genes):  
  relGenes = {}  
  for (gene, seg) in genes:
    if gene.info['data'][-4] not in gene_map:
      continue
    #print i.data
    tstart = gene.start
    tend = gene.end
    if not (gene.info['data'][-4].startswith("LOC") or gene.info['data'][-4].startswith("LINC")):
      if gene not in relGenes:
        relGenes[gene] = [tstart, tend, [(int(z[0]),int(z[1])) for z in zip(gene.info['data'][9].split(','), gene.info['data'][10].split(',')) if z[0] != '']]
      else:
        oldTStart = relGenes[gene][0]
        oldTEnd = relGenes[gene][1]
        if tstart < oldTStart:
          oldTStart = tstart
        if tend > oldTEnd:
          oldTEnd = tend        
        relGenes[gene] = (oldTStart,oldTEnd,relGenes[2]+[(int(z[0]),int(z[1])) for z in zip(gene.info['data'][9].split(','), gene.info['data'][10].split(',')) if z[0] != ''])
  return relGenes

def parse_genes(chrom):
  print("Building interval tree for chrom " + chrom)
  t = IntervalTree()
  with open("/pedigree2/projects/namphuon/programs/docker/aa/data_repo//hg19/human_hg19_september_2011/Genes_July_2010_hg19.gff") as infile:
    for line in infile:
      fields = line.rsplit("\t")
      #gene = fields[-2]
      currChrom = fields[2]
      tstart = int(fields[4])
      tend = int(fields[5])
      if chrom == currChrom:
        t[tstart:tend] = fields
  return t

#need to add functionality for plotting exon posns
def plot_gene_track(currStart, genes, segment, total_length_with_spacing, strand, exons = True):
  delta = control_delta if 'control' in segment.info else 0
  for gene in genes:    
    if gene.info['Name'].upper().startswith('RP11-') or gene.info['Name'].startswith('RP4-') or gene.info['Name'].startswith('RNU6') or gene.info['Name'].startswith('SNOR') or gene.info['Name'].startswith('U6') or gene.info['Name'].startswith('AC099') or gene.info['Name'].startswith('RP5') or gene.info['Name'].startswith('U4') or gene.info['Name'].find('RNA') != -1:
      continue
    #e_posns is a list of tuples of exon (start,end)
    #do exons
    tstart,tend = (gene.start,gene.end)
    #tstart,tend = (ensembl_map[gene.info['data'][-4]].start,ensembl_map[gene.info['data'][-4]].end)
    seg_len = segment.end - segment.start
    if strand == "+":
      normStart = currStart - max(0,tstart-segment.start)
      normEnd = currStart - min(segment.end-segment.start,tend-segment.start)
    else:
      normEnd = currStart - min(segment.end-segment.start,segment.end-tstart)
      normStart = currStart - max(0,segment.end - tend)
    # if tstart < pTup[1]:
    #   truncStart = True
    # if tend > pTup[2]:
    #   truncEnd = True
    gstart_angle = normStart/total_length_with_spacing*360
    gend_angle = normEnd/total_length_with_spacing*360
    
    text_angle = (gstart_angle + gend_angle)/2.0
    if gend_angle < 0 and gstart_angle > 0:
      gend_angle+=360      
    
    patches.append(mpatches.Wedge((0,0), outer_bar+delta, gend_angle, gstart_angle, width=bar_width/2.0))    
    color = 'r' if gene.info['Name'].upper() in oncogene_map else 'teal'    
    f_color_v.append(color)
    e_color_v.append(color)
    lw_v.append(0)
    
    x,y = pol2cart(outer_bar+(bar_width/2.0),(text_angle/360*2*np.pi))
    x_t,y_t = pol2cart(outer_bar + bar_width + gene_delta + delta,(text_angle/360*2*np.pi))
    #ax.plot([x,x_t],[y,y_t],color='grey',linewidth=0.4)
    
    if text_angle < -90 and text_angle > -360:
      text_angle+=180
      ax.text(x_t,y_t,gene.info['Name'],color=color,rotation=text_angle,
        ha='right',fontsize=font_size,rotation_mode='anchor')    
    else:
       ax.text(x_t,y_t,gene.info['Name'],color=color,rotation=text_angle,
        ha='left',fontsize=font_size,rotation_mode='anchor')
    #ax.plot([x,x_t],[y,y_t],color='grey',linewidth=0.2)  

def plot_exon_track(currStart, genes, segment, total_length_with_spacing, strand, color='black'):
  #   gene_names = set([gene_map[gene.info['data'][-4]].upper() for gene in genes])
  #   exons = [e for e in hg19.exon_list if e.info['Name'] in gene_names]
  #   exon_names = set([e.info['Name'] for e in exons])
  #   missing = [gene for gene in genes if gene_map[gene.info['data'][-4]].upper() not in exon_names and not (gene_map[gene.info['data'][-4]].upper()[0:5] == 'RP11-' or gene_map[gene.info['data'][-4]].upper()[0:4] == 'RP4-' or gene_map[gene.info['data'][-4]].upper()[0:4] == 'RNU6')]
  #   for m in missing:
  #     starts = m.info['data'][9].split(',')[0:-1]
  #     ends = m.info['data'][10].split(',')[0:-1]
  #     for i in xrange(0,len(starts)):
  #       exons.append(hg19.interval(m.chrom, int(starts[i]),int(ends[i])))
  delta = control_delta if 'control' in segment.info else 0
  for gene in genes: 
    exons = [g for g in gene.info['intervals'] if 'exon_number' in g.info['data']]
    for e in exons:
      if e.intersection(segment) is None:
        continue
      tstart,tend = (e.start,e.end)
      seg_len = segment.end - segment.start
      if strand == "+":
        normStart = currStart - max(0,tstart-segment.start)
        normEnd = currStart - min(segment.end-segment.start,tend-segment.start)
      else:
        normEnd = currStart - min(segment.end-segment.start,segment.end-tstart)
        normStart = currStart - max(0,segment.end - tend)
      gstart_angle = normStart/total_length_with_spacing*360
      gend_angle = normEnd/total_length_with_spacing*360
    
      if gend_angle < 0 and gstart_angle > 0:
        gend_angle+=360
    
      patches.append(mpatches.Wedge((0,0), outer_bar+delta, gend_angle, gstart_angle, width=bar_width/0.85))    
      f_color_v.append(color)
      e_color_v.append(color)
      lw_v.append(0)

#g_locs,cycle[cycle_num],seg_padding
def plot_ref_genome(start_points,lens,cycle,total_length_with_spacing,seg_padding,bed_feat_list=[]):
    seg_posns = []
    previous_end = total_length_with_spacing*(global_rot/360.0)
    for ind,sp in enumerate(start_points):
      start_point = int(previous_end - sp)
      start_angle = start_point/total_length_with_spacing*360
      end_angle = (start_point - lens[ind])/total_length_with_spacing*360
    
      #segseqD referenced as global variable here because I'm lazy
      seg_coord_tup = segSeqD[cycle[ind][0]]
      text_angle = (start_angle + end_angle)/2.0
      # text_angle = start_angle
      delta = control_delta if 'control' in seg_coord_tup.info else 0
      x,y = pol2cart(outer_bar + bar_width*-5+delta,text_angle/360*2*np.pi)
    
      #Writes the really ugly label for the segment on the outside
      #segment_name = "Segment " + cycle[ind][0] + cycle[ind][1] + "\n" + seg_coord_tup.chrom + ":" + str(seg_coord_tup.start) + "-" + str(seg_coord_tup.end)
      #segment_name = "Segment " + cycle[ind][0] + cycle[ind][1] + "\n" + seg_coord_tup.chrom
      segment_name = seg_coord_tup.chrom
      if text_angle < -90 and text_angle > -360:
        text_angle-=180
        ha = "right"
      else:
        ha = "left"
      ax.text(x,y,segment_name,color='k',rotation = text_angle,ha=ha,fontsize=font_size-1,rotation_mode='anchor')
    
      #makes the reference genome wedges  
      if end_angle < 0 and start_angle > 0:
        end_angle+=360
      
      patches.append(mpatches.Wedge((0,0), outer_bar+delta, end_angle, start_angle, width=bar_width))
      #f_color_v.append(chromosome_colors[seg_coord_tup.chrom])
      f_color_v.append('silver' if 'color' not in seg_coord_tup.info else seg_coord_tup.info['color'])
      e_color_v.append('silver')
      lw_v.append(0.2)
    
      #makes the ticks on the reference genome wedges
      if cycle[ind][1] == "+":
        posns = zip(range(seg_coord_tup.start,seg_coord_tup.end+1),np.arange(start_point,start_point-lens[ind]-1,-1))
      else:
        posns = zip(np.arange(seg_coord_tup.end,seg_coord_tup.start-1,-1),np.arange(start_point,start_point-lens[ind]-1,-1))
      for j in posns:
        if j[0] % tick_positions == 0 and j[0] != seg_coord_tup.start and j[0] != seg_coord_tup.end:
          text_angle = j[1]/total_length_with_spacing*360
          x,y = pol2cart(outer_bar+delta,(text_angle/360*2*np.pi))
          xp,yp = pol2cart(outer_bar+delta+minor_tick_delta/2.,(text_angle/360*2*np.pi))
          x_t,y_t = pol2cart(outer_bar+delta+minor_tick_delta,(text_angle/360*2*np.pi))
          ax.plot([x,xp],[y,yp],color='grey',linewidth=0.2)
        
          if text_angle < -90 and text_angle > -360:
            text_angle-=180
            ha = "right"
            txt = str(int(round((j[0])/10000))) + " "
          else:
            ha = "left"
            txt = " " + str(int(round((j[0])/10000)))
          ax.text(x_t,y_t,txt,color='grey',rotation=text_angle,
          ha=ha,fontsize=font_size/1.5,rotation_mode='anchor')    
        if j[0] == seg_coord_tup.start or j[0] == seg_coord_tup.end:
          text_angle = j[1]/total_length_with_spacing*360
          x,y = pol2cart(outer_bar+delta,(text_angle/360*2*np.pi))
          xp,yp = pol2cart(outer_bar+delta+minor_tick_delta/2.,(text_angle/360*2*np.pi))
          x_t,y_t = pol2cart(outer_bar+delta-bar_width,(text_angle/360*2*np.pi))
          ax.plot([x,xp],[y,yp],color='grey',linewidth=0.2)  
          if text_angle < -90 and text_angle > -360:
            text_angle-=180
            ha = "right"
            txt = str(j[0]) + " "
          else:
            ha = "left"
            txt = " " + str(j[0])
          ax.text(x_t,y_t,txt,color='black',rotation=text_angle,
          ha=ha,fontsize=font_size/2+1,rotation_mode='anchor')    
      #overhaul
      genes = [h[0] for h in ensembl_grc37.intersection([seg_coord_tup])]
      genes = [g for g in genes if (not g.info['Name'].upper().startswith('RP11-') and not g.info['Name'].startswith('RP4-') and not g.info['Name'].startswith('RNU6') or g.info['Name'].startswith('SNOR') and not g.info['Name'].startswith('U6') or g.info['Name'].startswith('AC099') and not g.info['Name'].startswith('RP5'))]
      plot_gene_track(start_point,genes,seg_coord_tup,total_length_with_spacing,cycle[ind][1])
      plot_exon_track(start_point,genes,seg_coord_tup,total_length_with_spacing,cycle[ind][1])

#OVERHAUL
def plot_ref_cmaps(start_point,seg_cmap_vector):
  start_angle = start_point/total_length_with_spacing*360
  end_angle = (start_point - seg_cmap_vector[-1])/total_length_with_spacing*360
  if end_angle < 0 and start_angle >0:
    end_angle+=360  
  patches.append(mpatches.Wedge((0,0), mid_bar+bar_width, end_angle, start_angle, width=bar_width))
  f_color_v.append('darkorange')
  e_color_v.append('k')
  lw_v.append(0)
  
  lab_locs = []
  print seg_cmap_vector
  for i in seg_cmap_vector[:-1]:
    label_angle = (start_point - i)/total_length_with_spacing*2*np.pi
    x,y = pol2cart(mid_bar,label_angle)
    x_t,y_t = pol2cart(mid_bar+bar_width,label_angle)
    ax.plot([x,x_t],[y,y_t],color='k',alpha=0.9,linewidth=0.2)
    lab_locs.append(start_point - i)
  return lab_locs

#parse the alignment file
def parse_alnfile(path_aln_file):
  aln_vect = []
  with open(path_aln_file) as infile:
    meta_header = infile.next().rstrip()[1:].split()
    aln_metadata_fields = infile.next().rstrip()[1:].split()
    meta_dict = dict(zip(meta_header,aln_metadata_fields))
    aln_header = infile.next().rstrip()[1:].split()
    for line in infile:
      fields = line.rstrip().split()
      fields_dict = dict(zip(aln_header,fields))
      aln_vect.append(fields_dict)
  return aln_vect,meta_dict

#parse cycles file
def parse_cycles_file(cycles_file, addchr = False):
  cycles = {}
  segSeqD = {}
  with open(cycles_file) as infile:
    for line in infile:
      if line.startswith("Segment"):
        fields = line.rstrip().split()
        lowerBound = int(fields[3])
        upperBound = int(fields[4])
        chrom = fields[2]
        if addchr:
          chrom = "chr%s" % chrom
        segNum = fields[1]
        segSeqD[segNum] = hg19.interval(chrom,lowerBound,upperBound, info={'name':segNum})
      elif "Cycle=" in line:
        curr_cycle = []
        fields = line.rstrip().rsplit(";")
        lineD = {x.rsplit("=")[0]:x.rsplit("=")[1] for x in fields}
        segs = lineD["Segments"].rsplit(",")
        #TODO: Need to rotate the segs in case there's a 0 within the path, or else path is incorrect
        skip = False
        for i in segs:
          seg = i[:-1]
          if seg != "0":
            strand = i[-1]
            curr_cycle.append((seg,strand))
          else:
            skip = True
        if not skip:
          cycles[lineD["Cycle"]] = curr_cycle
  return cycles,segSeqD

#parse bed file
def parse_bed_file(bed_file):
  bed_list = []
  print "f",bed_file
  with open(bed_file) as infile:
    for line in infile:
      fields = line.rstrip().split()
      fields[1] = int(fields[1])
      fields[2] = int(fields[2])
      fields[3] = float(fields[3])
      bed_list.append(fields)
  return bed_list

def get_contig_locs(aln_vect,meta_dict):
  #go through the seg_seq and calculate total lengths of alignments
  seg_seq = meta_dict["seg_seq"].split(",")
  total_seg_length = [seg_cmap_lens[i] for i in seg_seq]
  contig_list = []
  prev = None
  for a_d in aln_vect:
    if a_d["contig_id"] != prev:
      prev = a_d["contig_id"]
      contig_list.append[prev]
  total_contig_length = sum([contig_cmap_lens[x] for i in contig_set])
  total_genomic_length = max(total_seg_length,total_contig_length)
  spacing = contig_spacing*total_contig_length
  total_genomic_length+=(spacing*len(contig_set))
  
  start_points = [0]
  for ind in len(contig_list)-1:
    start_points.append(start_points[-1] + spacing + contig_cmap_lens[contig_list[ind]]) #this gets the previous ind!
  return dict(zip(contig_list,start_points))

def get_seg_locs_from_contig_aln():
  pass

#get start locations for a cycle
def get_seg_locs_from_cycle(cycle,segSeqD):
  lens = []
  for i in cycle:
    curr_len = segSeqD[i[0]].end - segSeqD[i[0]].start
    lens.append(curr_len)
  total_seg_len = sum(lens)
  seg_padding = seg_spacing*total_seg_len
  start_points = [0]
  for i in lens[:-1]:
    start_points.append(start_points[-1] + i + seg_padding)
  total_length = start_points[-1] + lens[-1] + seg_padding
  return start_points,lens,float(total_length), seg_padding

def feat_bed_to_lookup(bed_list):
  bed_dict = defaultdict(list)
  for i in bed_list:
    if i[0].startswith("hs"):
      i[0] = i[0].replace("hs","chr")
    bed_dict[i[0]].append(i[1:])
  return bed_dict

def find_peaks(bed_data, window = 1000):
  bed_data.sort()
  current = bed_data[0]
  keeps = hg19.interval_list()
  start = current
  for i in xrange(0,len(bed_data)):
    if bed_data[i].chrom != start.chrom:
      keeps.append(current)
      current = bed_data[i]
      start = current      
    elif bed_data[i].intersects(start,10000) and bed_data[i].info['value'] > current.info['value']:
      current = bed_data[i]
    elif not bed_data[i].intersects(start,window):
      keeps.append(current)
      current = bed_data[i]
      start = current
  keeps.append(current)
  return keeps

# if __name__ == '__main__':
# Parses the command line arguments
parser = argparse.ArgumentParser(description="Corrects and extends alignment paths to produce BioNano contig/AA segment scaffolds")
parser.add_argument("--bionano_alignments",help="turns on Bionano visualization mode (requires contigs,segs,key,path_alignment args")
parser.add_argument("-c", "--contigs", help="contig cmap file")
parser.add_argument("-s", "--segs", help="segments cmap file")
parser.add_argument("--cycles_file",help="cycles file")
parser.add_argument("--cycle",help="cycle number to visualize")
parser.add_argument("-k", "--keyfile", help="segments cmap key file")
parser.add_argument("-i", "--path_alignment", help="AR path alignment file")
parser.add_argument("--sname", help="output prefix")
parser.add_argument("--bed_files",help="bed file list",nargs="+")
parser.add_argument("--feature_labels",help="bed feature names",nargs="+")
args = parser.parse_args()
args.cycles_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/COLO320_DM_S270/onco_amplicon1_cycles.txt'
args.cycle = '16'
args.fpkm_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/COLO320_DM_S270/colo320dm.fpkm.csv'

if not args.sname:
  args.sname = args.segs.rsplit(".")[0] + "_"
else:
  samp_name = args.sname.rsplit("/")[-1]

fname = samp_name

bed_feat_dict = {}
if args.bed_files:
  for i,j in zip(args.bed_files,args.feature_labels):
    print j,i
    #feature name -> chromosome -> ordered list of positions
    bed_list = parse_bed_file(i)
    bed_feat_dict[j] = feat_bed_to_lookup(bed_list)

outer_bar = max(bed_track_height*(len(bed_feat_dict)+2),10)

bed_data = hg19.interval_list([hg19.interval('chr8', 127638302, 127938302, info={'value':int(random.random()*100)}), hg19.interval('chr8', 128716346,128746346, info={'value':int(random.random()*100)})])
bed_data.sort()

args.prefix_name = '/pedigree2/projects/namphuon/programs/CycleViz/COLO320DM'
args.cycles_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/COLO320_DM_S270/onco_amplicon1_cycles.txt'
args.fpkm_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/COLO320_DM_S270/colo320dm.fpkm.csv'
args.wgs_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/COLO320_DM_S270/colo320dm.wgs.1000.pileup.log.bed'
cycles_numbers = ['6', '9', '10', '12', '13', '14', '15', '16','19']
args.atac_peak_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/ATAC-seq/SRC1655_summits_250ext_q1e6_nochrM_merged.bed'
args.atac_file = '/pedigree2/projects/namphuon/results/paul_gbm39/ATAC/COLO320DM.atac.1000.pileup.log.bed'


args.prefix_name = '/pedigree2/projects/namphuon/programs/CycleViz/PC3'
args.cycles_file = '/nucleus/pedigree/projects/extrachromosome/data/turner2017/reconstruction/run14/FF-77_amplicon4_cycles.txt'
args.fpkm_file = '/pedigree2/projects/namphuon/results/paul_gbm39/rnaseq/PC3.fpkm.csv'
args.wgs_file = '/pedigree2/projects/namphuon/results/paul_gbm39/PC3/PC3.wgs.1000.pileup.log.bed'
args.atac_file = '/pedigree2/projects/namphuon/results/paul_gbm39/ATAC/PC3.atac.1000.pileup.log.bed'
args.atac_peak_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/ATAC-seq/YD_320_summits_250ext_q1e6_nochrM_merged.bed'
cycles_numbers = ['3', '9','10', '13', '32']

args.prefix_name = '/pedigree2/projects/namphuon/programs/CycleViz/GBM39KT'
args.cycles_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/GBM39-KT_S274/onco_amplicon1_cycles.txt'
args.fpkm_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/GBM39-KT_S274/GBM39-KT_S274.fpkm.csv'
args.wgs_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/GBM39-KT_S274/GBM39-KT_S274.wgs.1000.pileup.log.bed'
args.atac_file = '/pedigree2/projects/namphuon/results/paul_gbm39/ATAC/GBM39KT.atac.1000.pileup.log.bed'
args.atac_peak_file = '/pedigree2/projects/namphuon/data/paul_gbm39/unsorted/ATAC-seq/YD_316_summits_250ext_q1e6_nochrM_merged.bed'
args.fourc_file = '/pedigree2/projects/namphuon/results/paul_gbm39/4C/ramya/TL317_final_norm.bedgraph'
args.fourc_base_file = '/pedigree2/projects/namphuon/results/paul_gbm39/4C/ramya/TL125_final_norm.bedgraph'
primer = hg19.interval('chr7',55088102,55088151)
cycles_numbers = ['4','261']

args.prefix_name = '/pedigree2/projects/namphuon/programs/CycleViz/TCGA-L7-A6VZ'
args.cycles_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-L7-A6VZ-01A-12_WGS133091339321_CNV21180_amplicon5_cycles.txt'
args.fpkm_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-L7-A6VZ.fpkm.bed'
args.wgs_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-L7-A6VZ.wgs.1000.pileup.log.bed'
cycles_numbers = ['34']
args.atac_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-L7-A6VZ.atac.bed'
cycle_num = '34'

args.prefix_name = '/pedigree2/projects/namphuon/programs/CycleViz/TCGA-A7-A0D9'
args.cycles_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-A7-A0D9-01A-31_WGS24514959190_CNV20421_amplicon1_cycles.txt'
args.fpkm_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-A7-A0D9.fpkm.bed'
args.wgs_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-A7-A0D9.wgs.1000.pileup.log.bed'
args.atac_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-A7-A0D9.atac.bed'
cycles_numbers = [8]
cycle_num = '8'

args.prefix_name = '/pedigree2/projects/namphuon/programs/CycleViz/TCGA-AO-A124'
args.cycles_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-AO-A124-01A-11_WGS143473425566_CNV49327_amplicon14_cycles.txt'
args.fpkm_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-AO-A124.fpkm.bed'
args.wgs_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-AO-A124.wgs.1000.pileup.log.bed'
args.atac_file = '/pedigree2/projects/namphuon/results/paul_gbm39/howard/TCGA-AO-A124.atac.bed'
cycles_numbers = [1]
cycle_num = '1'

for cycle_num in cycles_numbers:  
  cycles,segSeqD = parse_cycles_file(args.cycles_file,addchr=True)
  foo = [(cid, sum([segSeqD[c[0]].end-segSeqD[c[0]].start for c in cycle])) for (cid, cycle) in cycles.items()]
  foo = sorted(foo, key=lambda x: -x[1])
  
  fname = '%s.%s' % (args.prefix_name, cycle_num)
  plt.clf()
  fig, ax = plt.subplots()
  
  patches = []
  f_color_v = []
  e_color_v = []
  lw_v = []
  cycle = cycles[cycle_num]
  segment_cycles = hg19.interval_list([segSeqD[c[0]] for c in cycle])
  segment_cycles.sort()
  #if cycle_num == '12' and args.prefix_name =='/pedigree2/projects/namphuon/programs/CycleViz/COLO320DM':
  #  cycle.insert(1,('start','+'))
  #  #cycle.insert(1,('end','+'))    
  #  segSeqD['start'] = hg19.interval('chr8',127019489,127039489,info={'color':'pink','control':True})
  #  #segSeqD['end'] = hg19.interval('chr8',129003683,129023683,info={'color':'pink'})      
  #if cycle_num == '4' and args.prefix_name =='/pedigree2/projects/namphuon/programs/CycleViz/GBM39KT':
  #  cycle.insert(2,('start','+'))
  #  cycle.insert(2,('end','+'))    
  #  segSeqD['start'] = hg19.interval('chr7',54810974,54830975,info={'color':'pink','control':True})
  #  segSeqD['end'] = hg19.interval('chr7',56117061,56137061,info={'color':'pink','control':True})    
  #if cycle_num == '3' and args.prefix_name =='/pedigree2/projects/namphuon/programs/CycleViz/PC3':
  #  cycle.insert(4,('start','+'))
  #  #cycle.insert(4,('end','+'))
  #  segSeqD['start'] = hg19.interval('chr8',112231503,112251503,info={'color':'pink','control':True})
  #  #segSeqD['end'] = hg19.interval('chr8',129545706,129565705,info={'color':'pink'})
  #if cycle_num == '34' and args.prefix_name =='/pedigree2/projects/namphuon/programs/CycleViz/TCGA-L7-A6VZ':
  #  cycle.insert(2,('start','+'))
  #  #cycle.insert(4,('end','+'))
  #  segSeqD['start'] = hg19.interval('chr12',24220000,24300000,info={'color':'pink','control':True})
  #if cycle_num == '8' and args.prefix_name =='/pedigree2/projects/namphuon/programs/CycleViz/TCGA-A7-A0D9':
  #  cycle.insert(2,('start','+'))
  #  #cycle.insert(4,('end','+'))
  #  segSeqD['start'] = hg19.interval('chr20',63000020,63025520,info={'color':'pink','control':True})
  start_points,lens,total_length_with_spacing,seg_padding = get_seg_locs_from_cycle(cycle,segSeqD)
  plot_ref_genome(start_points, lens, cycle,total_length_with_spacing,seg_padding,args.feature_labels)
  
  p = PatchCollection(patches)
  p.set_facecolor(f_color_v)
  p.set_edgecolor(e_color_v)
  p.set_linewidth(lw_v)
  ax.add_collection(p)
  ax.set_xlim(-(outer_bar+5), (outer_bar+5))
  ax.set_ylim(-(outer_bar+5), (outer_bar+5))
  ax.set_aspect(1.0)
  
  #wgs_bed = BuildBedTrack(0, BuildBedTrack.load_bed(args.wgs_file,log=True), is_log = False, num_axis=5,ymin=0,ymax=8)
  beds = BuildBedTrack.load_bed(args.wgs_file,log=True)  
  beds.sort()  
  for b in beds:
    b.info['value']*=2
  max([b.info['value'] for b in beds])
  foo = [segSeqD[s[0]] for s in cycle]  
  foo.sort()  
  foo = beds.intersection(foo)
  wgs_bed = BuildBedTrack(0, beds, is_log = False, num_axis=2, color='orange',ymin = 0, ymax = 125, track_pos = (17.5,14))
  wgs_bed.build_axis(line_one = 2)
  axis = PatchCollection(wgs_bed.bed_patches,match_original=True)
  ax.add_collection(axis)
  wgs_bed.build_segments()
  
  fpkm = load_fpkm_file(args.fpkm_file)
  fpkm = dict([(k,float(v)) for (k,v) in fpkm.items() if k in ensembl_grc37_map])
  fpkm_bed = build_genebed_from_fpkm(fpkm)
  #for f in fpkm_bed:
  #  f.info['value'] = 10**f.info['value']
  rna_bed = BuildBedTrack(1, fpkm_bed, is_log = True, ymin = 1e-3, ymax = 1000, num_axis=5, color='red', track_pos = (12.5,9))
  #rna_bed.build_axis(line_one = 1,ticks=[1e-3,1e-2,1e-1,1,10,100,1000])
  rna_bed.build_axis(line_one = 1,ticks=[1e-3,1e-2,1e-1,1,10,100,1000])
  axis = PatchCollection(rna_bed.bed_patches,match_original=True)
  ax.add_collection(axis)
  rna_bed.build_segments()
  
  #plt.axis('off')
  #plt.savefig(fname + '.png',dpi=600)
  #plt.close()
  #print "done"
  
  #   atac_peak = BuildBedTrack.load_bed(args.atac_peak_file, value = 1)
  atac_bed =  BuildBedTrack.load_bed(args.atac_file,log=False)
  foo = find_peaks(atac_bed,window=10000)
  foo.sort()
  
  #atac_bed = hg19.interval_list([h for h in atac_bed if h.info['value'] >= 10])
  #   for a in atac_peak:
  #     a.info['color'] = 'red'
  atac_seq = BuildBedTrack(2, foo, is_log = True, ymin = 1e-3, ymax = 1000, num_axis=5, color='green',point_spacing = 1000,track_pos = (6,-2))
  atac_seq.track_rmax=atac_seq.track_rmax+1
  #hits = [h[0] for h in atac_seq.bed_data.intersection(atac_peak)]
  #for h in hits:
  #  h.info['fill'] = 'solid'
  #  h.info['color'] = 'green'  
  atac_seq.build_axis(line_one = True,ticks=[1,10,100,1000])
  #atac_seq.bed_patches = hg19.interval_list(hits)
  axis = PatchCollection(atac_seq.bed_patches,match_original=True)
  ax.add_collection(axis)
  atac_seq.build_segments()
  #atac_seq.plot_highlight(start_points, cycles[cycle_num], bed_data = atac_peak)
  
  plt.axis('off')
  plt.savefig(fname + '.pdf', bbox_inches='tight')


class BuildBedTrack:
  
  def __init__(self, track, bed_data, ymax = None, ymin = None, is_log = False, num_axis = 10, point_spacing = 500, color = 'k', color_bed = None, track_pos = None):
    self.track = track
    if track_pos is None:
      self.track_rmin = outer_bar-bar_width-((track+1)*bed_track_height)+bed_track_spacing-1
      self.track_rmax = outer_bar-bar_width-((track)*bed_track_height)-bed_track_spacing-1
    else:
      (self.track_rmax, self.track_rmin) = (track_pos[0],track_pos[1])
    self.is_log = is_log
    self.ymax = max([s.info['value'] for s in bed_data]) if ymax is None else ymax
    self.ymin = min([s.info['value'] for s in bed_data]) if ymin is None else ymin
    if self.is_log and self.ymin == 0:
      self.ymin = 1e-3
    self.num_axis = num_axis
    self.point_spacing = point_spacing
    self.bed_patches = []
    self.color = color
    self.bed_data = bed_data
    self.color_bed = color_bed
  
  def build_axis(self, line_one = None, ticks = None):
    if ticks is None:
      i = self.num_axis*1.    
      while (i >= 0):
        y_value = self.ymin+(self.ymax-self.ymin)/(self.num_axis*1.)*i
        if self.is_log:
          y_value = math.log10(self.ymin)+(math.log10(self.ymax)-math.log10(self.ymin))/(self.num_axis*1.)*i
          y_value = math.log10(y_value)
        r_value = self.track_rmin+(self.track_rmax-self.track_rmin)/(self.num_axis)*i
        self.bed_patches.append(mpatches.Wedge((0, 0), r_value, 0, 360, width=0.01))
        x_t,y_t = pol2cart(r_value,np.pi*2*1/4)
        #ax.text(x_t,y_t,int(y_value) if self.is_log == False else "%.0E" % Decimal(y_value),color='grey', ha='right',fontsize=3,rotation_mode='anchor')
        ax.text(x_t,y_t,int(y_value) if self.is_log == False else "%0.1f" % Decimal(10**y_value),color='grey', ha='right',fontsize=font_size/1.5,rotation_mode='anchor')
        i+=-1  
    else:
      for t in ticks:
        if self.is_log:
          y_value = (math.log10(t)-math.log10(self.ymin))/(math.log10(self.ymax)-math.log10(self.ymin))
        else:
          y_value = self.ymin+(self.ymax-self.ymin)/(self.num_axis*1.)*t
        r_value = y_value*(self.track_rmax-self.track_rmin)+self.track_rmin
        self.bed_patches.append(mpatches.Wedge((0, 0), r_value, 0, 360, width=0.01))
        x_t,y_t = pol2cart(r_value,np.pi*2*1/4)
        ax.text(x_t,y_t,"%0.3f" % t,color='grey', ha='right',fontsize=font_size/1.5,rotation_mode='anchor')        
    if line_one is not None:
      if self.is_log:
        y_scale_value = (math.log10(line_one)-math.log10(self.ymin))/ (math.log10(self.ymax)-math.log10(self.ymin))
      else:
        y_scale_value = (line_one-self.ymin)/(self.ymax-self.ymin)          
      r_scale_value = y_scale_value*(self.track_rmax-self.track_rmin)+self.track_rmin
      self.bed_patches.append(mpatches.Wedge((0, 0), r_scale_value, 0, 360, alpha=0.5, width=0.025, edgecolor = None, facecolor = 'darkblue', color = 'darkblue', fill=False, clip_on= False))             
  
  def draw_beziar_curve(self, primer, bed_data = None, fudge = 0, lw = 0.25):
    if bed_data is None:
      bed_data = self.bed_data
    points_x = []
    points_y = []
    strength = []
    color = []
    previous_end = total_length_with_spacing*(global_rot/360.0)
    start_x = None
    start_y = None
    for ind,sp in enumerate(start_points):
      start_point = int(previous_end - sp)
      start_angle = start_point/total_length_with_spacing*360
      end_angle = (start_point - lens[ind])/total_length_with_spacing*360
      
      #segseqD referenced as global variable here because I'm lazy
      segment = segSeqD[cycle[ind][0]] 
      strand = cycle[ind][1]
      hits = [h[0] for h in bed_data.intersection([segment])]
      if start_x is None and primer.intersection(segment) is not None:
        if strand == "+":
          normStart = start_point - max(0,primer.start-segment.start)
          normEnd = start_point - min(segment.end-segment.start,primer.start-segment.start)
        else:
          normEnd = start_point - min(segment.end-segment.start,segment.end-primer.start)
          normStart = start_point - max(0,segment.end - primer.start)
        start_x,start_y = pol2cart(self.track_rmin+self.delta,normStart/total_length_with_spacing*2*np.pi)        
      for h in hits:
        pos = int((h.start+h.end)/2)
        if pos > segment.end or pos < segment.start:
          continue
        if strand == "+":
          normStart = start_point - max(0,pos-segment.start)
          normEnd = start_point - min(segment.end-segment.start,pos-segment.start)
        else:
          normEnd = start_point - min(segment.end-segment.start,segment.end-pos)
          normStart = start_point - max(0,segment.end - pos)
        x_s,y_s = pol2cart(self.track_rmin+self.delta,(normStart+fudge)/total_length_with_spacing*2*np.pi)
        foo = points_x.append(x_s)
        foo = points_y.append(y_s)
        foo = strength.append(h.info['value'])
        foo = color.append(h.info['color'])
    for x in xrange(0,len(points_x)):
      plt.plot([start_x,points_x[x]], [start_y,points_y[x]], linewidth=lw, color=color[x],alpha=0.5)
      #pp1 = mpatches.PathPatch(
      #Path([(start_x, start_y), (0, 0), (points_x[x], points_y[x]), (start_x, start_y)],
      #     [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CLOSEPOLY]),
      #fc="none", transform=ax.transData)
      #ax.add_patch(pp1)
    
  def build_segments(self, bed_data = None):
    if bed_data is None:
      bed_data = self.bed_data
    points_x = []
    points_y = []
    colors = []
    
    fpoints_x = []
    fpoints_y = []
    fcolors = []
    
    previous_end = total_length_with_spacing*(global_rot/360.0)    
    for ind,sp in enumerate(start_points):
      start_point = int(previous_end - sp)
      start_angle = start_point/total_length_with_spacing*360
      end_angle = (start_point - lens[ind])/total_length_with_spacing*360
      
      #segseqD referenced as global variable here because I'm lazy
      segment = segSeqD[cycle[ind][0]] 
      strand = cycle[ind][1]
      hits = [h[0] for h in bed_data.intersection([segment])]
      if self.color_bed is not None:
        color_subhits = hg19.interval_list([h[0] for h in self.color_bed.intersection([segment])])      
      for h in hits:
        for pos in xrange(h.start, h.end, self.point_spacing):
          if pos > segment.end or pos < segment.start:
            continue
          if self.color_bed is not None:
            temp = hg19.interval(h.chrom, pos, pos)
            color_hits = color_subhits.intersection([temp],self.point_spacing)
            if len(color_hits) != 0:
              color = color_hits[0][0].info['color']
            else:
              color = self.color if 'color' not in h.info else h.info['color']
          else:
            color = self.color if 'color' not in h.info else h.info['color']
          if strand == "+":
            normStart = start_point - max(0,pos-segment.start)
            normEnd = start_point - min(segment.end-segment.start,pos-segment.start)
          else:
            normEnd = start_point - min(segment.end-segment.start,segment.end-pos)
            normStart = start_point - max(0,segment.end - pos)
          hvalue = h.info['value'] if h.info['value'] > self.ymin else self.ymin
          hvalue = hvalue if hvalue < self.ymax else self.ymax
          y_scale_value = (1.*hvalue-self.ymin)/(self.ymax-self.ymin)          
          if self.is_log:
            y_scale_value = (math.log10(hvalue)-math.log10(self.ymin))/(math.log10(self.ymax)-math.log10(self.ymin))
          r_scale_value = y_scale_value*(self.track_rmax-self.track_rmin)+self.track_rmin                    
          x_s,y_s = pol2cart(r_scale_value,normStart/total_length_with_spacing*2*np.pi)
          if 'fill' in h.info:
            foo = fpoints_x.append(x_s)
            foo = fpoints_y.append(y_s)
            fcolors.append(color)          
          else:
            foo = points_x.append(x_s)
            foo = points_y.append(y_s)
            colors.append(color)
    foo = ax.scatter(points_x,points_y,marker='o',s=1,linewidths=0.01,facecolors='none',color=colors)  
    foo = ax.scatter(fpoints_x,fpoints_y,marker='*',s=1,linewidths=0.01,color=fcolors)  
  
  @staticmethod
  def load_bed(bed_file, value = None, log = False, sep='\t'):
    bed_data = hg19.interval_list()
    for line in open(bed_file):
      res = line.split(sep)
      if value is None:
        bed_data.append(hg19.interval(res[0], int(res[1]), int(res[2]), info={'value':float(res[3]) if not log else 10**float(res[3])}))
      else:
        bed_data.append(hg19.interval(res[0], int(res[1]), int(res[2]), info={'value':value}))
    bed_data.sort()
    return bed_data
  
  def plot_highlight(self, start_points, cycle, bed_data = None):
    if bed_data is None:
      bed_data = self.bed_data
    seg_posns = []
    previous_end = total_length_with_spacing*(global_rot/360.0)    
    for ind,sp in enumerate(start_points):
      segment = segSeqD[cycle[ind][0]]    
      strand = cycle[ind][1]
      start_point = int(previous_end - sp)
      beds = [h[0] for h in bed_data.intersection([segment])]
      for bed in beds:    
        tstart,tend = bed.start, bed.end
        seg_len = segment.end - segment.start
        if strand == "+":
          normStart = start_point - max(0,tstart-segment.start)
          normEnd = start_point - min(segment.end-segment.start,tend-segment.start)
        else:
          normEnd = start_point - min(segment.end-segment.start,segment.end-tstart)
          normStart = start_point - max(0,segment.end - tend)
        gstart_angle = normStart/total_length_with_spacing*360
        gend_angle = normEnd/total_length_with_spacing*360    
        if gend_angle < 0 and gstart_angle > 0:
          gend_angle+=360
        patches.append(mpatches.Wedge((0,0), self.track_rmax, gend_angle, gstart_angle, width=self.track_rmax-self.track_rmin, alpha=0.6))    
        color = 'gray'
        f_color_v.append(color)
        e_color_v.append(color)
        lw_v.append(0)

if False:
  cycles,segSeqD = parse_cycles_file(args.cycles_file,addchr=False)
  foo = [(cid, sum([segSeqD[c[0]].end-segSeqD[c[0]].start for c in cycle])) for (cid, cycle) in cycles.items()]
  foo = sorted(foo, key=lambda x: -x[1])
  cycle_num = '4'
  fname = '%s.%s' % (args.prefix_name, cycle_num)
  plt.clf()
  fig, ax = plt.subplots()
  
  patches = []
  f_color_v = []
  e_color_v = []
  lw_v = []
  cycle = cycles[cycle_num]  
  #cycle.insert(2,('start','+'))
  #cycle.insert(2,('end','+'))    
  #segSeqD['start'] = hg19.interval('chr7',54810974,54830975,info={'color':'pink','control':True})
  #segSeqD['end'] = hg19.interval('chr7',56117061,56137061,info={'color':'pink','control':True})    
  
  start_points,lens,total_length_with_spacing,seg_padding = get_seg_locs_from_cycle(cycle,segSeqD)
  plot_ref_genome(start_points, lens, cycle,total_length_with_spacing,seg_padding,args.feature_labels)
  
  p = PatchCollection(patches)
  p.set_facecolor(f_color_v)
  p.set_edgecolor(e_color_v)
  p.set_linewidth(lw_v)
  ax.add_collection(p)
  ax.set_xlim(-(outer_bar+5), (outer_bar+5))
  ax.set_ylim(-(outer_bar+5), (outer_bar+5))
  ax.set_aspect(1.0)
  
  #   beds = BuildBedTrack.load_bed(args.wgs_file,log=True)  
  #   beds.sort()
  #   foo = [segSeqD[s[0]] for s in cycle]
  #   foo.sort()
  #   foo = beds.intersection(foo)
  #   wgs_bed = BuildBedTrack(0, beds, is_log = False, num_axis=5, color='blue',ymin = 0, ymax = 20)
  #   wgs_bed.build_axis()
  #   axis = PatchCollection(wgs_bed.bed_patches)
  #   ax.add_collection(axis)
  #   wgs_bed.build_segments()
  
  base_4c = BuildBedTrack.load_bed(args.fourc_base_file, sep = ' ')
  normal_4c = BuildBedTrack.load_bed(args.fourc_file, sep=' ') 
  normals = ['/pedigree2/projects/namphuon/results/paul_gbm39/4C/ramya/TL318_final_norm.bedgraph','/pedigree2/projects/namphuon/results/paul_gbm39/4C/ramya/TL319_final_norm.bedgraph']
  for n in normals:
    ns = BuildBedTrack.load_bed(n, sep=' ') 
    for i in xrange(0,len(base_4c)):
      if normal_4c[i].start != ns[i].start:
        print "failed"
        break
      normal_4c[i].info['value']+=ns[i].info['value']
  for n in normal_4c:
    n.info['value']/=3              
  remove = [] 
  for i in xrange(0,len(base_4c)):
    if normal_4c[i].intersection(base_4c[i]) is None:
      print "Broke"
      break
    else:
      if normal_4c[i].info['value'] == 0 and base_4c[i].info['value'] == 0:
        remove.append(normal_4c[i])
        continue
      normal_4c[i].info['value'] = (normal_4c[i].info['value'])/(base_4c[i].info['value']+normal_4c[i].info['value'])
  normal_4c = hg19.interval_list([v for v in normal_4c if v not in remove])
  
  keeps = ['/pedigree2/projects/namphuon/results/paul_gbm39/4C/ramya/TL125_nearbait_highinter.bed', '/pedigree2/projects/namphuon/results/paul_gbm39/4C/ramya/TL317_nearbait_highinter.bed', '/pedigree2/projects/namphuon/results/paul_gbm39/4C/ramya/TL318_nearbait_highinter.bed','/pedigree2/projects/namphuon/results/paul_gbm39/4C/ramya/TL319_nearbait_highinter.bed']
  keepers = hg19.interval_list()
  for k in keeps:
    ns = BuildBedTrack.load_bed(k, value=1) 
    ns.sort()
    keepers.extend(ns)
  keepers.sort()
  normal_4c.sort()
  normal_4c = hg19.interval_list(set([h[0] for h in normal_4c.intersection(keepers)]))
  gradient = np.linspace(0, 1, 256)
  colors = plt.cm.get_cmap('Reds',100)
  x = np.linspace(0.0, 1.0, 101)
  edge_map = colors(x)[np.newaxis, :, :3]
  for n in normal_4c:
    n.info['color'] = edge_map[0][max(int(100*n.info['value']**2),1)]  
    #n.info['color'] = edge_map[0][int((n.info['value'])*100)+0]  
  four_seq = BuildBedTrack(0, normal_4c, is_log = True, ymin = 1e-3, ymax = 100000, num_axis=5, color='blue')
  four_seq.track_rmin=outer_bar-2
  #four_seq.build_axis(ticks=[1,10,100,1000,10000,100000])
  axis = PatchCollection(four_seq.bed_patches)
  ax.add_collection(axis)
  #four_seq.build_segments(base_4c)    
  four_seq.delta = 1
  four_seq.color = 'blue'   
  #four_seq.draw_beziar_curve(primer, bed_data = base_4c)    
  four_seq.color = 'red'        
  #four_seq.build_segments(normal_4c)  
  four_seq.draw_beziar_curve(primer, bed_data = normal_4c, lw=1)   
  plt.axis('off')  
    
  #cmap = matplotlib.cm.Reds
  #norm = matplotlib.colors.Normalize(vmin=0, vmax=100)
  #cb1 = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap,
  #                                norm=norm,
  #                                orientation='horizontal')
  #cb1.set_label('4C ratio')  
  plt.savefig(fname +'.4c'+ '.pdf', bbox_inches='tight')
  #plt.savefig(fname + '.png',dpi=600)
  plt.close()
  print "done"
  
  fname = '%s.%s.color' % (args.prefix_name, cycle_num)
  # Make plot with vertical (default) colorbar
  fig, ax = plt.subplots()
  data = np.clip(randn(250, 250), 0, 1)
  cax = ax.imshow(data, interpolation='nearest', cmap=cm.Reds)
  ax.set_title('Gaussian noise with vertical colorbar')
  # Add colorbar, make sure to specify tick locations to match desired ticklabels
  cbar = fig.colorbar(cax, ticks=[-1, 0, 1])
  cbar.ax.set_yticklabels(['< -1', '0', '> 1'])  # vertically oriented colorbar
  # Make plot with horizontal colorbar
  fig, ax = plt.subplots()
  data = np.clip(randn(250, 250), 0, 1)
  cax = ax.imshow(data, interpolation='nearest', cmap=cm.Reds)
  ax.set_title('Gaussian noise with horizontal colorbar')
  cbar = fig.colorbar(cax, ticks=[0, .25**2, .5**2, .75**2,1], orientation='vertical')
  cbar.ax.set_yticklabels(['0.0', '0.25', '0.5', '0.75','1.0'],fontdict={'fontsize':12})  # horizontal colorbar
  plt.show()
  plt.show()
  plt.savefig(fname +'.4c'+ '.pdf', bbox_inches='tight')
  #plt.savefig(fname + '.png',dpi=600)
  plt.close()
  print "done"

