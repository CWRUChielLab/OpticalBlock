from neuron import h
import pylab # (used for plotting)

# A NEURON model of an unmylenated axon
class Axon:
    # NEURON defines a coordinate system along the length of a segment from 0 to 1;
    # for consistency we will define 0 as the left side of the axon and 1 as the
    # right side.
    left_side = 0.
    right_side = 1.
    middle = 0.5

    # create the axon
    def __init__(self,
            num_sections = 1000, # number of compartments along the axon's length
            length = 5000.,      # total length of axon (um)
            diam = 1.,           # axonal diameter (um)
            ):

        # store some parameters for future use
        self.length = length
        self.stimuli = []

        # create the sections
        self.sections = [h.Section() for i in range(num_sections)]

        # initialize the sections
        for sec in self.sections:
            # set the geometry
            sec.L = length/num_sections
            sec.diam = diam

            # set the passive properties
            sec.cm = 1.
            sec.Ra = 100.

            # insert the active channels
            #sec.insert('hh')
            sec.insert('fhm1')

            sec.localtemp_fhm1 = 16

        # connect the sections end to end into a chain, connecting the right
        # side of each segment to the left side of the next
        for i in range(num_sections - 1):
            self.sections[i].connect(self.sections[i+1], Axon.left_side,
                    Axon.right_side)

    # get the section a given fraction of the length along the axon
    def section_at_f(self, f):
        if (f == 1):
            # corner case: normally a section extends from the left boundary of
            # the section up to, but *not* including, the right boundary (which
            # is the left boundary of the next section).  The rightmost section
            # needs to include its right boundary, however.
            return self.sections[-1]
        else:
            return self.sections[int(f * len(self.sections))]

    # get the section containing the given position along the axon
    def section_at_x(self, x):
        return self.section_at_f(x / self.length)

    # insert a simple current at the given position
    def insert_stim(self, x=0.):
        stim = h.IClamp(Axon.middle, self.section_at_x(x))
        stim.amp = 1.  # nA
        stim.delay = 1. # ms
        stim.dur = 1.   # ms

        # NOTE: NEURON will remove the stimulus as soon as it's no longer
        # reachable from python code, so we need to store it.
        self.stimuli.append(stim)


# if we're running this code directly (vs. importing it as a library),
# run a simple simulation and generate a demo plot
if __name__ == "__main__":
    axon = Axon()
    axon.insert_stim()

    # set up recording vectors
    t = h.Vector()
    t.record(h._ref_t)

    num_v_traces = 10
    v_traces = []
    for i in range(num_v_traces):
        v = h.Vector()
        v.record(axon.section_at_f(i/(num_v_traces - 1.))
            (Axon.middle)._ref_v)
        v_traces.append(v)

    # initialize the simulation
    h.dt = 0.025 # integration time step, in ms
    tstop = 19.9   # duration of integration
    v_init = -65 # initial membrane potential, in mV
    h.finitialize(v_init)
    h.fcurrent()

    # run the simulation
    while h.t < tstop:
        h.fadvance()

    # plot the results
    pylab.figure(1, figsize=(6,6))
    for v in v_traces:
        pylab.plot(t, v)
    pylab.xlabel("time (ms)")
    pylab.ylabel("membrane potential (mV)")
    pylab.show()

    # pdf and png export seem to be broken in NEURON - see comment here
    # http://www.neuron.yale.edu/phpbb/viewtopic.php?f=2&t=2097 -
    # thus we export the figure as an eps
    pylab.savefig("demo_plot.eps") #, bbox_inches='tight')
