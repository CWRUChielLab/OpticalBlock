from neuron import h
import pylab # (used for plotting)

# A NEURON model of an unmylenated axon
class Axon:
    # NEURON defines a coordinate system along the length of a section from 0 to 1;
    # for consistency we will define 0 as the left side of the axon and 1 as the
    # right side.
    left_side = 0.
    right_side = 1.
    middle = 0.5


    # Create the axon
    def __init__(self,
            num_sections = 1000, # number of compartments along the axon's length
            length = 1000.,      # total length of axon (um)
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


    # Get the index of the section at the given length along the axon
    def section_id_at_f(self, f):
        # parameter check
        if f < 0.:
            raise ValueError("position beyond left end of axon")
        if f > 1.:
            raise ValueError("position beyond right end of axon")

        if (f == 1):
            # corner case: normally a section extends from the left boundary of
            # the section up to, but *not* including, the right boundary (which
            # is the left boundary of the next section).  The rightmost section
            # needs to include its right boundary, however.
            return len(self.sections) - 1
        else:
            return int(f * len(self.sections));


    # Get the section a given fraction of the length along the axon
    def section_at_f(self, f):
        return self.sections[self.section_id_at_f(f)]


    # Get the section containing the given position along the axon
    def section_at_x(self, x):
        return self.section_at_f(x / self.length)


    # Apply a function to all the sections that are part of a given range
    # of positions along the axon.  The end sections of the range are
    # included, even though they may be only partially within the range.
    def apply_to_sections(self,
            func,         # function to apply, of the form func(sec, x_center)
            x_start = 0., # start of the range of positions, measured from left
            x_end = None  # end of the range; if None everything to the right
                          # x_start will be included.
            ):
        # parameter checking and cleanup
        if x_end == None:
            x_end = self.length
        if x_end < x_start:
            raise ValueError("x_end must be to the right of x_start")

        for i in range(self.section_id_at_f(x_start / self.length),
                self.section_id_at_f(x_end / self.length) + 1):
            func(self.sections[i],
                    (i + Axon.middle) * self.length / len(self.sections))


    # Insert a simple current at the given position
    def insert_stim(self, x=0.):
        stim = h.IClamp(Axon.middle, self.section_at_x(x))
        stim.amp = 1.    # nA
        stim.delay = 0.1 # ms
        stim.dur = 1.    # ms

        # NOTE: NEURON will remove the stimulus as soon as it's no longer
        # reachable from python code, so we need to store it.
        self.stimuli.append(stim)


    # Set the temperature of a range of sections
    def set_temp(self,
            temp,       # a temperature (Celsius) or a function of the form
                        # g(x) that returns the temperature at position x.
            x_start=0,  # start of the range of positions (0 == left end)
            x_end=None  # end of the range of positions (None == right end)
            ):
        # if temperature is not a function, turn it into one.
        if hasattr(temp, '__call__'):
            temp_at_x = temp
        else:
            temp_at_x = lambda x: temp

        def set_section_temp(sec, x):
            sec.localtemp_fhm1 = temp_at_x(x)

        self.apply_to_sections(set_section_temp, x_start, x_end);



# if we're running this code directly (vs. importing it as a library),
# run a simple simulation and generate a demo plot
if __name__ == "__main__":
    axon = Axon()
    axon.insert_stim()
    axon.set_temp(54.49, 300, 600)

    # set up recording vectors
    t = h.Vector()
    t.record(h._ref_t)

    num_v_traces = 20
    v_traces = []
    for i in range(num_v_traces):
        v = h.Vector()
        v.record(axon.section_at_f(
            # record at num_v_traces points along the axon, equally spaced
            # from eachother and from the end points (since we don't care
            # about things like the impedance mismatch at the ends)
            (i + 1) * 1.0 / (num_v_traces + 1))
            (Axon.middle)._ref_v)
        v_traces.append(v)

    # initialize the simulation
    h.dt = 0.005 # integration time step, in ms
    tstop = 2.99 # duration of integration
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
