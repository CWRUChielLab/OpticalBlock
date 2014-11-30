from neuron import h
import pylab # (used for plotting)

# NEURON defines a coordinate system along the length of a segment from 0 to 1;
# for consistency we will define 0 as the left side of the axon and 1 as the
# right side.
left_side = 0.
right_side = 1.

# create the axon
def create_axon(
        num_sections = 1000, # number of compartments along the axon's length
        length = 1000.,      # total length of axon (um)
        diam = 1.,           # axonal diameter (um)
        ):

    # create the sections
    sections = [h.Section() for i in range(num_sections)]

    # initialize the sections
    for sec in sections:
        # set the geometry
        sec.L = length/num_sections
        sec.diam = diam

        # set the passive properties
        sec.cm = 1.
        sec.Ra = 100.

        # insert the active channels
        sec.insert('hh')

    # connect the sections end to end into a chain, connecting the right
    # side of each segment to the left side of the next
    for i in range(num_sections - 1):
        sections[i].connect(sections[i+1], left_side, right_side)

    return sections

axon = create_axon()

# insert the stimulus
stim = h.IClamp(left_side, axon[0])
stim.amp = 10.
stim.delay = 5.
stim.dur = 5.

# set up recording vectors
t = h.Vector()
t.record(h._ref_t)
v = h.Vector()
v.record(axon[-1](right_side)._ref_v)

# run the simulation
h.load_file("stdrun.hoc")
h.init()
h.tstop = 20.
h.run()

# plot the results
pylab.figure(1, figsize=(6,6))
pylab.plot(t, v)
pylab.xlabel("time (ms)")
pylab.ylabel("membrane potential (mV)")
pylab.show()

# pdf and png export seem to be broken in NEURON - see comment here
# http://www.neuron.yale.edu/phpbb/viewtopic.php?f=2&t=2097 -
# thus we export the figure as an eps
pylab.savefig("plot.eps") #, bbox_inches='tight')
