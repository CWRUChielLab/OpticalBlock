from neuron import h

import re   # regular expressions
import json # used for reading the config file
import sys  # used for command line parsing
import bisect
import copy
import math

# A NEURON model of an unmylenated axon
class Axon:
    # NEURON defines a coordinate system along the length of a section from 0 to 1;
    # for consistency we will define 0 as the left side of the axon and 1 as the
    # right side.
    left_side = 0.
    right_side = 1.
    middle = 0.5


    # Create the axon
    def __init__(self, config, delay_config=False):

        # store parameters for future use
        self.config = config
        self.stimuli = []

        # create the sections
        self.sections = [h.Section() for i in range(config['num_sections'])]

        # connect the sections end to end into a chain, connecting the right
        # side of each segment to the left side of the next
        for i in range(len(self.sections) - 1):
            self.sections[i].connect(self.sections[i+1], Axon.left_side,
                    Axon.right_side)

        # insert the active channels
        for sec in self.sections:
            #sec.insert('hh')
            sec.insert('fhm1')

        if not delay_config:
            # set the various section parameters
            self.update_sections()


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
        return self.section_at_f(x / float(config['axon_length']))


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
            x_end = float(config['axon_length'])
        if x_end < x_start:
            raise ValueError("x_end must be to the right of x_start")

        for i in range(self.section_id_at_f(x_start / float(config['axon_length'])),
                self.section_id_at_f(x_end / float(config['axon_length'])) + 1):
            func(self.sections[i],
                    (i + Axon.middle) * float(config['axon_length']) / len(self.sections))


    # Insert a simple current at the given position
    def insert_stim(self, x=0.):
        stim = h.IClamp(Axon.middle, self.section_at_x(x))
        stim.amp = 1.    # nA
        stim.delay = 0.1 # ms
        stim.dur = 1.    # ms

        # NOTE: NEURON will remove the stimulus as soon as it's no longer
        # reachable from python code, so we need to store it.
        self.stimuli.append(stim)


    def set_section_temp(self, sec, temp):
        sec.localtemp_fhm1 = temp
        sec.Ra = (self.config['axial_resistance'] *
                self.config['axial_resistance_Q10']**(
                (temp - self.config['axial_resistance_T'])/10.))


    # updates the section parameters based on the current config
    def update_sections(self):
        def set_section_params(sec, x):
            # fill in the position and evaluate any configuration settings
            # that use it.
            config = self.config.copy()
            config[u'x'] = x
            config = simplify_config(config)

            # set the geometry
            sec.L = float(config['axon_length'])/len(self.sections)
            sec.diam = config['axon_diameter']

            # set the passive properties
            sec.cm = config['membrane_capacitance']

            # set the active properties
            sec.pnabar_fhm1 = config['g_Na_bar']
            sec.ppbar_fhm1 = config['g_Na_persistent_bar']
            sec.pkbar_fhm1 = config['g_K_bar']

            # set the temperature
            #print x, config['axon_temperature']
            self.set_section_temp(sec, config['axon_temperature'])

        self.apply_to_sections(set_section_params);


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
            self.set_section_temp(sec, temp_at_x(x))

        self.apply_to_sections(set_section_temp, x_start, x_end);


    # Set the sodium conductance of a range of sections
    def set_gNa(self,
            gNa,        # a conductance (S/cm^2( or a function of the form
                        # g(x) that returns the conductance at position x.
            x_start=0,  # start of the range of positions (0 == left end)
            x_end=None  # end of the range of positions (None == right end)
            ):
        # if conductance is not a function, turn it into one.
        if hasattr(gNa, '__call__'):
            gNa_at_x = gNa
        else:
            gNa_at_x = lambda x: gNa

        def set_section_gNa(sec, x):
            sec.pnabar_fhm1 = gNa_at_x(x)

        self.apply_to_sections(set_section_gNa, x_start, x_end);



# perform a binary search to bound where a function switches from true to
# false at least once (e.g from not blocked to blocked as a function of
# temperature).  Returns a pair containing the new upper and lower bounds.
# Throws a ValueError if the function has the same value at the upper and
# lower bounds.
def boolean_bisect(
        func,          # a function of the form f(x) that returns a boolean
        lower_bound,   # a value known to be below the change
        upper_bound,   # a value known to be above the change
        num_iterations # the number of bisections to do to narrow the bounds
        ):
    xlow = lower_bound
    xhigh = upper_bound

    upperval = func(xhigh)
    lowerval = func(xlow)
    if (upperval == lowerval):
        raise ValueError(
                "Function has the same value at the upper and lower bounds")

    for i in range(num_iterations):
        xmid = (xlow + xhigh) / 2.
        if func(xmid) == upperval:
            xhigh = xmid
        else:
            xlow = xmid

    return (xlow, xhigh)


# perform linear interpolation given the list of points describing the
# function (xs and ys) and a point at which to evaluate the function, x
def interpolate(xs, ys, x):
    # assume the function is constant beyond the end points given
    if math.isnan(x):
        return float('nan')
    elif x <= xs[0]:
        return ys[0]
    elif x >= xs[-1]:
        return ys[-1]
    else:
        i = bisect.bisect(xs, x)
        f = float(xs[i] - x)/(xs[i] - xs[i-1])
        return ys[i-1] * f + ys[i] * (1 - f)


def is_numeric(x):
    return type(x) == int or type(x) == float

def is_numeric_list(x):
    if type(x) != list:
        return False
    for item in x:
        if not is_numeric(item):
            return False
    return True


def simplify_config(config):

    newconfig = copy.deepcopy(config)

    def simplify_pass(value, context):
        if type(value) == str or type(value) == unicode:
            # replace variables with their valueues
            if (context.has_key(value)):
                return context[value], True
        elif type(value) == list:
            list_changed = False
            for i in range(len(value)):
                value[i],item_changed = simplify_pass(value[i], context)
                list_changed = list_changed or item_changed
            return value, list_changed
        elif type(value) == dict:
            dict_changed = False

            # simplify the components
            for key,val in value.items():
                value[key], item_changed = simplify_pass(value[key], context)
                dict_changed = dict_changed or item_changed

            # see if we can collapse the whole thing
            if value.has_key('action'):
                action = value['action']
                if action == 'interpolate':
                    if (is_numeric_list(value['example_inputs']) and
                            is_numeric_list(value['example_outputs']) and
                            is_numeric(value['new_input'])):
                        return (
                            interpolate(value['example_inputs'],
                                value['example_outputs'], value['new_input']),
                            True)
                elif action == 'interpolate_from_csv':
                    if is_numeric(value['new_input']):
                        # read the csv file
                        vals = []
                        with open(value['csv_file'], 'r') as f:
                            for line in f.readlines():
                                vals.append([float(val) for val in line.split(',')])
                        example_inputs, example_outputs = zip(*vals)
                        return (interpolate(example_inputs, example_outputs,
                            value['new_input']),
                            True)
                elif action == 'gaussian':
                    if (is_numeric(value['center']) and
                            is_numeric(value['width']) and
                            is_numeric(value['height']) and
                            is_numeric(value['input'])):
                        # calculate the gaussian
                        result = value['height'] * math.exp(
                            -float(value['input'] - value['center'])**2 /
                            value['width']**2)
                        return (result, True)


            return value, dict_changed

        return value, False

    changed = True
    while changed:
        newconfig,changed = simplify_pass(newconfig, newconfig)

    return newconfig


# runs the model and returns True iff an action potential reaches the end
# (last 1%) of the axon
def is_blocked(axon):
    v = h.Vector()
    v.record(axon.section_at_f(0.99)
        (Axon.middle)._ref_v)

    # initialize the simulation
    h.dt = 0.005 # integration time step, in ms
    tstop = 2.99 # duration of integration
    v_init = -65 # initial membrane potential, in mV
    h.finitialize(v_init)
    h.fcurrent()

    # run the simulation
    while h.t < tstop:
        h.fadvance()

    threshold = -45 # mV
    if max(v) >= threshold:
        return False
    else:
        return True



def record_plot(
        func,
        xmin,
        xmax,
        title,
        xlabel,
        ylabel,
        num_points=11
        ):
    csv_filename = title.replace(' ','_') + '.csv'

    xs = []
    ys = []
    for i in range(num_points):
        # space the points equally from max_width to min_width
        x = i * (xmax - xmin) * 1.0 / (num_points - 1) + xmin
        try:
            y = func(x)
        except ValueError:
            y = float("NaN")
        print((x, y))
        xs.append(x)
        ys.append(y)

    # save the data as a csv
    with open(csv_filename, 'w') as csv_file:
        # write a header
        csv_file.write("{0},{1}\n".format(
                xlabel.replace(' ','_').replace('(','').replace(')',''),
                ylabel.replace(' ','_').replace('(','').replace(')','')
                ))

        # write each data point
        for x, y in zip(xs, ys):
            csv_file.write("{0},{1}\n".format(x, y))


def run_single_simulation(config, interactive):
    axon = Axon(config)
    axon.insert_stim()

    # set up recording vectors for python plots and the csv file
    t = h.Vector()
    t.record(h._ref_t)

    num_v_traces = 30
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

    # set up NEURON plotting code (if we're in an interactive session)
    if interactive:
        g = h.Graph()
        g.size(0, 3, -80, 55)
        for i in range(num_v_traces):
            g.addvar('v(0.5)',
                    sec=axon.section_at_f((i+1) * 1.0 / (num_v_traces + 1)))

    # initialize the simulation
    h.dt = config['max_time_step']
    tstop = config['integration_time']
    h.finitialize(config['initial_membrane_potential'])
    h.fcurrent()

    # run the simulation
    if interactive:
        g.begin()
        while h.t < tstop:
            h.fadvance()
            g.plot(h.t)
        g.flush()
    else:
        while h.t < tstop:
            h.fadvance()

    # save the data as a csv
    with open('demo_traces.csv', 'w') as csv_file:
        # start with a header of the form "t_ms, V0_mV, V1_mv, V2_mV,..."
        csv_file.write(", ".join(
            ["t_ms"] + ["V{0}_mV".format(i) for i in range(num_v_traces)]
            ) + "\n")

        # write the time and each of the recorded voltages at that time
        for row in zip(t, *v_traces):
            csv_file.write(", ".join([str(x) for x in row]) + "\n")

# searches a config setting to see if it depends on the given variable
def has_variable(expr, var):
    if type(expr) == str or type(expr) == unicode:
        return expr == var
    elif type(expr) == list or type(expr) == tuple:
        for item in expr:
            if has_variable(item, var):
                return True
        return False
    elif type(expr) == dict:
        for key,val in expr.items():
            if has_variable(val, var):
                return True
        return False

def run_sweep_simulation(config, interactive):
    axon = Axon(config, delay_config=True)
    axon.insert_stim()

    # save the data as a csv
    csv_filename = "sweep.csv"
    with open(csv_filename, 'w') as csv_file:

        # Find the variables which change from sweep to sweep
        swept_vars = [key for key,val in config.items() if
                (has_variable(val, "sweep_param") or
                has_variable(val, "threshold_param")) and
                not has_variable(val, "x")]

        # write out the headers
        csv_file.write(", ".join(swept_vars) + "\n")

        for i in range(config['param_sweep_steps']):
            # space the points equally from max_width to min_width
            sweepconfig = copy.copy(config)
            sweepconfig["sweep_param"] = i * 1.0 / (
                    sweepconfig["param_sweep_steps"] - 1)
            sweepconfig = simplify_config(sweepconfig)

            def threshold_block_test(threshold_param):
                axon.config = sweepconfig.copy()
                axon.config[u"threshold_param"] = threshold_param
                axon.config = simplify_config(axon.config)
                axon.update_sections()
                print("  Testing " + ", ".join(
                    ["{0}:{1}".format(s,axon.config[s])
                        for s in swept_vars]) +
                    "...")
                blocked = is_blocked(axon)
                if blocked:
                    print("    blocked")
                else:
                    print("    not blocked")
                return blocked

            threshold_config = {}
            try:
                bounds = boolean_bisect(threshold_block_test, 0, 1, 10)
                threshold = sum(bounds)/2.

            except ValueError:
                # No threshold found
                threshold = float("NaN")

            # calculate all parameters at the threshold
            sweepconfig[u"threshold_param"] = threshold
            threshold_config = simplify_config(sweepconfig)

            # calculate all parameters at the threshold

            print("Threshold values: " + ", ".join(
                ["{0}:{1}".format(s,threshold_config[s])
                    for s in swept_vars]))

            # write out the values
            csv_file.write(", ".join(
                ["{0}".format(threshold_config[s])
                    for s in swept_vars]) + "\n")


# if we're running this code directly (vs. importing it as a library),
# run a simple simulation and generate a demo plot
if __name__ == '__main__':

    # We need to to some hacking to support command line parameters.  NEURON
    # tries to run all of the arguments on the command line, but we want to
    # use these extra parameters for config files.  This is not a disaster
    # because NEURON will try to run the config files after running all of the
    # code in this file, so errors don't prevent us from getting useful work
    # done (and JSON may not even cause an error).  They can, however, cause
    # NEURON to leave the user at a NEURON prompt, which prevents using this
    # code from makefiles, bash scripts, and other non-interactive tools.  To
    # work around this, if we're not running in interactive mode we can
    # terminate NEURON after running this script (before it tries to run the
    # remaining command line arguments).

    # Assume that the user doesn't want a NEURON prompt unless they've invoked
    # nrngui or specified '-' on the command line (standard for NEURON).
    interactive = "nrngui" in sys.argv[0] or "-" in sys.argv

    # extract the config files
    configfiles = (['defaultconfig.yaml'] +
        [arg for arg in sys.argv if ".yaml" in arg])

    print('Using config files: {0}\n'.format(", ".join(configfiles)))

    # read the configuration files in order, allowing later config files
    # to override earlier settings.
    config = {}
    for filename in configfiles:
        with open(filename,'r') as f:
            # read the entire file as text
            config_text = f.read(-1)

            # remove comments (from '#' to the end of the line)
            bare_config_text = re.sub('#[^\n]*\n', '\n', config_text)

            # parse it as JSON
            config.update(json.loads(bare_config_text))
    config = simplify_config(config)

    if config['param_sweep_steps'] == 1:
        run_single_simulation(config, interactive)
    else:
        run_sweep_simulation(config, interactive)

    # Now that we're done, quit if we're not in interactive mode
    # (see command line options above)
    if not interactive:
        h.quit();
