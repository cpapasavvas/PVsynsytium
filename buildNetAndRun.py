# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 23:40:57 2016

@author: osboxes
"""

import numpy
import simrun
import random
from neuron import h, gui
from neuronpy.util import spiketrain
from math import sin, cos, pi
from matplotlib import pyplot

class Cell(object):
    """Generic cell template."""
    def __init__(self):
        self.x, self.y, self.z = 0, 0, 0
        self.synlist = [] #### NEW CONSTRUCT IN THIS WORKSHEET
        self.create_sections()
        self.define_geometry()
        self.build_topology()
        self.build_subsets()
        self.define_biophysics()
        self.create_synapses()
    #
    def create_sections(self):
        """Create the sections of the cell. Remember to do this
        in the form::
            h.Section(name='soma', cell=self)
        """
        raise NotImplementedError("create_sections() is not implemented.")
    #
    def build_topology(self):
        """Connect the sections of the cell to build a tree."""
        raise NotImplementedError("build_topology() is not implemented.")
    #
    def define_geometry(self):
        """Set the 3D geometry of the cell."""
        raise NotImplementedError("define_geometry() is not implemented.")
    #
    def define_biophysics(self):
        """Assign the membrane properties across the cell."""
        raise NotImplementedError("define_biophysics() is not implemented.")
    #
    def create_synapses(self):
        """Subclasses should create synapses (such as ExpSyn) at various
        segments and add them to self.synlist."""
        pass # Ignore if child does not implement.
    #
    def build_subsets(self):
        """Build subset lists. This defines 'all', but subclasses may
        want to define others. If overridden, call super() to include 'all'."""
        self.all = h.SectionList()
        self.all.wholetree(sec=self.soma)
    #
    def connect2target(self, target, thresh=10):
        """Make a new NetCon with this cell's membrane
        potential at the soma as the source (i.e. the spike detector)
        onto the target passed in (i.e. a synapse on a cell).
        Subclasses may override with other spike detectors."""
        nc = h.NetCon(self.soma(0.5)._ref_v, target, sec = self.soma)
        nc.threshold = thresh
        return nc
    #
    def spikeDetector(self, thresh=0):
        """Make a new NetCon with this cell's membrane
        potential at the soma as the source (i.e. the spike detector).
        Subclasses may override with other spike detectors."""
        nc = h.NetCon(self.soma(0.5)._ref_v, None, sec = self.soma)
        nc.threshold = thresh
        return nc
    #
    def is_art(self):
        """Flag to check if we are an integrate-and-fire artificial cell."""
        return 0
    #
    def set_position(self, x, y, z):
        """
        Set the base location in 3D and move all other
        parts of the cell relative to that location.
        """
        for sec in self.all:
            for i in range(int(h.n3d())):
                h.pt3dchange(i,
                        x - self.x + h.x3d(i),
                        y - self.y + h.y3d(i),
                        z - self.z + h.z3d(i),
                        h.diam3d(i))
        self.x, self.y, self.z = x, y, z
    #
    def rotateZ(self, theta):
        """Rotate the cell about the Z axis."""
        rot_m = numpy.array([[sin(theta), cos(theta)], 
                              [cos(theta), -sin(theta)]])
        for sec in self.all:
            for i in range(int(h.n3d())):
                xy = numpy.dot([h.x3d(i), h.y3d(i)], rot_m)
                h.pt3dchange(i, xy[0], xy[1], h.z3d(i), h.diam3d(i))
                

class BallAndStick(Cell):  #### Inherits from Cell
    """Two-section cell: A soma with active channels and
    a dendrite with passive properties."""
    #### __init__ is gone and handled in Cell.
    #### We can override __init__ completely, or do some of
    #### our own initialization first, and then let Cell do its
    #### thing, and then do a bit more ourselves with "super".
    ####
    #### def __init__(self):
    ####     # Do some stuff
    ####     super(Cell, self).__init__()
    ####     # Do some more stuff
    #
    def create_sections(self):
        """Create the sections of the cell."""
        self.soma = h.Section(name='soma', cell=self)
        self.dendL = h.Section(name='dendL', cell=self)
        self.dendR = h.Section(name='dendR', cell=self)
        # self.axon = h.Section(name='axon', cell=self)
    #
    def build_topology(self):
        """Connect the sections of the cell to build a tree."""
        h.celsius = 23.0    # from Konstandoudaki (non physiological)
        self.dendL.connect(self.soma(1))
        self.dendR.connect(self.soma(1))
        # self.axon.connect(self.soma(0))
    #
    def define_geometry(self):
        """Set the 3D geometry of the cell."""
        self.soma.L = 27
        self.soma.diam = 29  # microns
        self.soma.nseg = 1
        # self.axon.L = 115
        # self.axon.diam = 1.5
        # self.axon.nseg = 1
        self.dendL.L = self.dendR.L = 200.0       # microns
        self.dendL.diam = self.dendR.diam = 0.8   # microns (Fukuda2003)
        self.dendL.nseg = self.dendR.nseg = 10
        # self.shape_3D()
    #
    def define_biophysics(self):
        """Assign the membrane properties across the cell."""
        #for sec in self.all: # 'all' exists in parent object.
        #    sec.Ra = 100    # Axial resistance in Ohm * cm
        #    sec.cm = 1      # Membrane capacitance in micro Farads / cm^2
        # Insert active Hodgkin-Huxley current in the soma
        #self.soma.insert('hh')
        #self.soma.gnabar_hh = 0.12  # Sodium conductance in S/cm2
        #self.soma.gkbar_hh = 0.036  # Potassium conductance in S/cm2
        #self.soma.gl_hh = 0.0003    # Leak conductance in S/cm2
        #self.soma.el_hh = -54.3     # Reversal potential in mV
        
        # notice that some values were originally different from the ones 
        # in the paper
        
        Rm = 10000
        
#        h.v_init = -70.7
#        my_ek = -89.4      # normal K state
        
        h.v_init = -48.9
        my_ek = -64.4       # high K state
        
        self.soma.insert('pas')
        self.soma.cm=1.2
        self.soma.g_pas=1.0 / Rm
#        self.soma.e_pas=-48.9
        self.soma.Ra=150.0
        
        self.soma.insert('Nafx')
        self.soma.gnafbar_Nafx= 0.045
        
        self.soma.insert('kdrin')           #delayed rectifier K+, S/cm2
        self.soma.gkdrbar_kdrin=0.018       #originally 0.018
        self.soma.ek = my_ek
        
        self.soma.insert('IKsin')           # D-type K+
        self.soma.gKsbar_IKsin=0.000725*0.1
        
        self.soma.insert('hin')
        self.soma.gbar_hin=0.00001
        
        self.soma.insert('kapin')            # A-type K+, S/cm2 
        self.soma.gkabar_kapin=0.0032*15     #originally 0.0032*15
        
        self.soma.insert('canin')
        self.soma.gcalbar_canin=0.0003
        
        self.soma.insert('kctin')           #fAHP
        self.soma.gkcbar_kctin=0.0001
        
        self.soma.insert('cadynin')
        
        # Insert passive current in the left dendrite
        self.dendL.insert('pas')
        self.dendL.cm=1.2
        self.dendL.g_pas=1.0 / Rm
        self.dendL.e_pas=-55.0
        self.dendL.Ra=150.0       
        
        self.dendL.insert('Nafx')
        self.dendL.gnafbar_Nafx=0.04
        
        self.dendL.insert('kdrin')          #delayed rectifier K+, S/cm2
        self.dendL.gkdrbar_kdrin=self.soma.gkdrbar_kdrin * 0.5
        self.dendL.ek = my_ek
        
        self.dendL.insert('kapin')          # A-type K+, S/cm2 
        self.dendL.gkabar_kapin= self.soma.gkabar_kapin*10
        
        # Insert passive current in the right dendrite
        self.dendR.insert('pas')
        self.dendR.cm=1.2
        self.dendR.g_pas=1.0 / Rm
        self.dendR.e_pas=-55.0
        self.dendR.Ra=150.0       
        
        self.dendR.insert('Nafx')
        self.dendR.gnafbar_Nafx=0.04
        
        self.dendR.insert('kdrin')          #delayed rectifier K+, S/cm2
        self.dendR.gkdrbar_kdrin=self.soma.gkdrbar_kdrin * 0.5
        self.dendR.ek = my_ek
        
        self.dendR.insert('kapin')          # A-type K+, S/cm2 
        self.dendR.gkabar_kapin= self.soma.gkabar_kapin*10
        
        # # Insert passive current in the axon
        # self.axon.insert('pas')
        # self.axon.cm = 1.2
        # self.axon.g_pas = 1.0 / Rm
        # self.axon.e_pas = -55.0
        # self.axon.Ra = 150.0
        
        # self.axon.insert('Nafx')
        # self.axon.gnafbar_Nafx= self.soma.gnafbar_Nafx * 10
        
        # self.axon.insert('kdrin')
        # self.axon.gkdrbar_kdrin = self.soma.gkdrbar_kdrin * 0.5
        
        
        self.current_balancein()
        

    #
    def shape_3D(self):
        """
        Set the default shape of the cell in 3D coordinates.
        Set soma(0) to the origin (0,0,0) and dend extending along
        the X-axis. -Innacurate
        """
        len1 = self.soma.L
        h.pt3dclear(sec=self.soma)
        h.pt3dadd(0, 0, 0, self.soma.diam, sec=self.soma)
        h.pt3dadd(0, len1, 0, self.soma.diam, sec=self.soma)
        len2 = self.dendL.L
        h.pt3dclear(sec=self.dendL)
        h.pt3dadd(0, len1, 0, self.dendL.diam, sec=self.dendL)
        h.pt3dadd(-len2*cos(pi/4),len2*sin(pi/4),0, 
                  self.dendL.diam,sec=self.dendL)
        h.pt3dclear(sec=self.dendR)
        h.pt3dadd(0, len1, 0, self.dendR.diam, sec=self.dendR)
        h.pt3dadd(len2*cos(pi/4),len2*sin(pi/4),0, 
                  self.dendR.diam,sec=self.dendR)
    #
    #### build_subsets, rotateZ, and set_location are gone. ####
    #
    #### NEW STUFF ####
    #
    #def create_synapses(self):
    #    """Add an exponentially decaying synapse in the middle
    #    of the soma. Set its tau to 2ms, and append this
    #    synapse to the synlist of the cell. It is used for injecting current"""
    #    syn = h.ExpSyn(self.soma(0.5))
    #    syn.tau = 2
    #    self.synlist.append(syn)
    #
    def current_balancein(self):
        """ This is a translation from Poirazi's current-balancein.hoc code
        """ 
        h.finitialize(-48.9)          # min 42
#        h.finitialize(-70.7)  
        h.fcurrent()
        
        for sec in self.all:
            if h.ismembrane('na_ion'):
                sec.e_pas = sec.v + sec.ina / sec.g_pas
            if h.ismembrane('k_ion'):
                sec.e_pas = sec.e_pas + sec.ik / sec.g_pas
            if h.ismembrane('ca_ion'):
                sec.e_pas = sec.e_pas + sec.ica / sec.g_pas
            if h.ismembrane('h'):
                sec.e_pas = sec.e_pas + sec.ihi / sec.g_pas
                
        h.fcurrent()
        
    
        
class Ring:
    """A network of *N* ball-and-stick cells where cell n makes an
    excitatory synapse onto cell n + 1 and the last, Nth cell in the
    network projects to the first cell.
    """
    def __init__(self, N=5, stim_w=0.04, stim_number=1,
            syn_w=0.01, syn_delay=5):
        """
        :param N: Number of cells.
        :param stim_w: Weight of the stimulus
        :param stim_number: Number of spikes in the stimulus
        :param syn_w: Synaptic weight
        :param syn_delay: Delay of the synapse
        """
        self._N = N              # Total number of cells in the net
        self.cells = []          # Cells in the net
        self.nclist = []         # NetCon list
        self.nclist_som = []     # NetCon list for spike detection
        self.stim = None         # Stimulator
        self.stim_w = stim_w     # Weight of stim
        self.stim_number = stim_number  # Number of stim spikes
        self.syn_w = syn_w       # Synaptic weight
        self.syn_delay = syn_delay  # Synaptic delay
        self.t_vec = h.Vector()   # Spike time of all cells
        self.id_vec = h.Vector()  # Ids of spike times
        self.set_numcells(N)  # Actually build the net.
    #
    def set_numcells(self, N, radius=50):
        """Create, layout, and connect N cells."""
        self._N = N
#        self.create_cells_chain(N)
        self.create_cells_slice(N)
        self.comp_distances()
        self.connectivityM()
#        self.connect_cells()   # obsolete
#        self.connect_stim_NetStim()  # obsolete
        self.connect_stim(200)
        self.spike_detector_list()
    #
    def create_cells_slice(self, N):
        """Create and layout N cells in the network."""
        self.cells = []
        #r = 50 # Radius of cell locations from origin (0,0,0) in microns
        N = self._N
        xlocations = 650 * numpy.random.uniform(0,1,N)
        xlocations = sorted(xlocations)
        
        for i in range(N):
            cell = BallAndStick()
            # When cells are created, the soma location is at (0,0,0) and
            # the dendrite extends along the X-axis.
            # First, at the origin, rotate about Z.
#            cell.rotateZ(-0.2)
            # Then reposition
            x_loc = xlocations[i]
            y_loc = 150 * numpy.random.uniform(0,1)
            z_loc = 150 * numpy.random.uniform(0,1)
            cell.set_position(x_loc, y_loc, z_loc)
            self.cells.append(cell)
            #
    def create_cells_chain(self, N):
        """Create and layout N cells in the network."""
        self.cells = []
        N = self._N
        for i in range(N):
            cell = BallAndStick()
            x_shift = 100
            cell.set_position(i*x_shift, 0, 0)
            self.cells.append(cell)
    #
    def connect_cells(self):
        """Connect cell n to cell n + 1. Obsolete"""
        self.nclist = []
        N = self._N
        for i in range(N):
            src = self.cells[i]
            tgt_syn = self.cells[(i+1)%N].synlist[0]
            nc = src.connect2target(tgt_syn)
            nc.weight[0] = self.syn_w
            nc.delay = self.syn_delay
            nc.record(self.t_vec, self.id_vec, i)
            self.nclist.append(nc)
    #
    def connect_stim_NetStim(self):
        """Connect a spiking generator to the first cell to get
        the network going. Obsolete"""
        self.stim = h.NetStim()
        self.stim.number = self.stim_number
        self.stim.start = 40
        self.ncstim = h.NetCon(self.stim, self.cells[0].synlist[0])
        self.ncstim.delay = 1
        self.ncstim.weight[0] = self.stim_w # NetCon weight is a vector.
    #
    def connect_stim(self, field):
        """Connect IClamp to all the cells up to x=field or only to the first 
         one if there is nothing up to x=field."""
        
        self.stimlist = []
        N = self._N
        for i in range(N):
            if self.cells[i].x < field:
                stim = h.IClamp(self.cells[i].soma(0.5))
                stim.dur = 25
                stim.amp = 0.6
                stim.delay = 60
                self.stimlist.append(stim)
            else:
                break
        
        if not self.stimlist:
            stim = h.IClamp(self.cells[0].soma(0.5))
            stim.dur = 25
            stim.amp = 0.6
            stim.delay = 60
            self.stimlist.append(stim)
            
    #
    def get_spikes(self):
        """Get the spikes as a list of lists."""
        return spiketrain.netconvecs_to_listoflists(self.t_vec, self.id_vec)
    #
    def conn_specif_cells(self, sec1, sec2, ind1, ind2, gap_res, flag):
        "electrically connects sec1 and sec2 reciprocally with weight"
        
        dist = self.distances[ind1][ind2]
        dendLength = sec1.L
        minpos = dist/dendLength/2.0
        if flag:
            pos = numpy.random.uniform(minpos,1)
        else:
            pos = minpos
        
        gapj1=h.gap(pos, sec=sec1)       
        gapj1.r=gap_res
        gapj2=h.gap(pos, sec=sec2)       
        gapj2.r=gap_res
        
        h.setpointer(sec2(pos)._ref_v, 'vgap', gapj1)
        h.setpointer(sec1(pos)._ref_v, 'vgap', gapj2)
        
        self.gaplist.append(gapj1)
        self.gaplist.append(gapj2)
#        self.nclist.append(nc1)
#        self.nclist.append(nc2)
        
    #
    def comp_distances(self):
        self.distances = [[0 for j in range(self._N)] for i in range(self._N)]
        
        for i in range(self._N-1):
            for j in range(i+1, self._N):
                x1 = self.cells[i].x
                x2 = self.cells[j].x
                y1 = self.cells[i].y
                y2 = self.cells[j].y
                z1 = self.cells[i].z
                z2 = self.cells[j].z
                self.distances[i][j]=((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5
    #
    def connectivityM(self):
        self.M = [[0 for j in range(self._N)] for i in range(self._N)]
        
        for i in range(self._N-1):
            for j in range(i+1, self._N):
                rndN = numpy.random.uniform(0,1)
                dist = self.distances[i][j]
                if rndN < (-0.6/400)*dist + 0.6:
                    self.M[i][j] = 1
                    self.M[j][i] = 1
    #
#    def connectivityM_Exp(self):
#        self.M = [[0 for j in range(self._N)] for i in range(self._N)]
#        
#        scaleN = 0.3
#        rndExp = numpy.random.exponential(4, self._N)   # do i really need this
#        rndExp = scaleN * (rndExp / max(rndExp))
#        for i in range(self._N-1):
#            maxPr = 1 - (scaleN - rndExp[i])
#            for j in range(i+1, self._N):
#                rndU = numpy.random.uniform(0,1)
#                dist = self.distances[i][j]
#                if rndU < (-1.0/400)*dist + maxPr:
#                    self.M[i][j] = 1
#                    self.M[j][i] = 1
    #
    def spike_detector_list(self):
        self.nclist_som = []
        N = self._N
        for i in range(N):
            src = self.cells[i]
            nc = src.spikeDetector()
            nc.record(self.t_vec, self.id_vec, i)
            self.nclist_som.append(nc)
        


myN= 70

for nSim in range(1000,1001):
        
    ring = Ring(N=myN)
    while ring.cells[20].x > 200.0:
        ring = Ring(N=myN)
        
    pyplot.imshow(ring.distances)
    pyplot.show()
    pyplot.imshow(ring.M)
    pyplot.show()
    pyplot.hist(map(sum, ring.M))
    
    gap_r = 1000       # this is given in MOhm (e.g., 5000 MOhm = 1/ (0.2 nS))
    
    ring.gaplist=[]
    for i in range(ring._N-1):
        for j in range(i+1, ring._N):
            if ring.M[i][j]:
                ring.conn_specif_cells(ring.cells[i].dendR, ring.cells[j].dendL,
                                       i, j, gap_r,1)
    
    ## this is for the chain
    #ring.gaplist=[]
    #for i in range(ring._N-1):
    #    ring.conn_specif_cells(ring.cells[i].dendR, ring.cells[i+1].dendL, 
    #                           i, i+1, gap_w, 0)
                
    
    soma_v_vec = [[] for y in range(myN)] 
    dend_v_vec = [[] for y in range(myN)] 
    t_vec = [[] for y in range(myN)] 
    for i in range(myN):
        soma_v_vec[i], dend_v_vec[i], t_vec[i] = \
        simrun.set_recording_vectors(ring.cells[i])
    
    simrun.simulate(tstop=100)
    
    for i in [5,15,25,35,45]:
        simrun.show_output(soma_v_vec[i], dend_v_vec[i], t_vec[i])
    pyplot.show()
    
    #Uncomment shape3D() above to make this worthwhile 
    #shape_window = h.PlotShape()
    #shape_window.exec_menu('Show Diam')
    
    spikes = ring.get_spikes()
    if len(spikes) > 1:
        from neuronpy.graphics import spikeplot
        filename = 'spikeplot'+str(nSim)+'.png'
        sp = spikeplot.SpikePlot(savefig=False)
        sp.set_fig_name(fig_name=filename)
        sp.set_marker('.')
        sp.set_markerscale(1.5)
        sp.set_markercolor('red')
        #sp.set_markeredgewidth(7)
        sp.plot_spikes(spikes)
