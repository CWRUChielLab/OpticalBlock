from neuron import h
import pylab # (used for plotting)

# create the axon
axon = h.Section()
axon.L = 1000.
axon.nseg = 100
axon.diam = 1.
axon.insert('hh')
axon.Ra = 100.
axon.cm = 1.
#for seg in axon:
#    seg.Ra = 100
#    seg.cm = 1

stim_position = 0.
stim = h.IClamp(stim_position, axon)
stim.amp = 10.
stim.delay = 5.
stim.dur = 5.

# set up recording vectors
t = h.Vector()
t.record(h._ref_t)
v = h.Vector()
v.record(axon(1.)._ref_v)

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
# thus we export the figure as a pdf
pylab.savefig("plot.eps") #, bbox_inches='tight')
