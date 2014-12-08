from nrnaxon import Axon
from neuron import h

# perform a binary search to bound where a function switches from true to
# false at least once (e.g from not blocked to blocked as a function of
# temperature).  Returns a pair containing the new upper and lower bounds.
# Throws a ValueError if the function has the same value at the upper and
# lower bounds.
def bisect(
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


def heat_block_temp_at_width(block_width):
    axon_length = 300.
    axon_num_sections = 3000
    axon = Axon(length=axon_length, num_sections=axon_num_sections)
    axon.insert_stim()

    def heat_block(temp):
        axon.set_temp(16)
        axon.set_temp(temp,
                (axon_length - block_width)/2.,
                (axon_length + block_width)/2.
                )

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

    bounds = bisect(heat_block, 0, 200, 20)
    return sum(bounds)/2.


if __name__ == "__main__":
    max_width = 82.
    min_width = 2.
    num_points = 20 + 1

    widths = []
    temps = []
    for i in range(num_points):
        # space the points equally from max_width to min_width
        width = i * (max_width - min_width) / (num_points - 1) + min_width
        try:
            temp = heat_block_temp_at_width(width)
        except ValueError:
            temp = float("NaN")
        print((width, temp))
        widths.append(width)
        temps.append(temp)

    # save the data as a csv
    with open('square_block_threshold.csv', 'w') as csv_file:
        # write a header
        csv_file.write("width_um,temp_C")

        # write each data point
        for width, temp in zip(widths, temps):
            csv_file.write("{0},{1}\n".format(width, temp))
