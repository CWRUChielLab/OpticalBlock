import nrnaxon
from neuron import h

def heat_block_temp_at_width(block_width):
    axon_length = 300.
    axon_num_sections = 3000
    axon = nrnaxon.Axon(length=axon_length, num_sections=axon_num_sections)
    axon.insert_stim()

    def heat_block(temp):
        axon.set_temp(16)
        axon.set_temp(temp,
                (axon_length - block_width)/2.,
                (axon_length + block_width)/2.
                )
        return nrnaxon.is_blocked(axon)

    bounds = nrnaxon.bisect(heat_block, 0, 200, 20)
    return sum(bounds)/2.


if __name__ == "__main__":
    nrnaxon.record_plot(
        func=heat_block_temp_at_width,
        xmin = 2,
        xmax = 82,
        title = "Block temperature vs square block width",
        xlabel = "Block width (um)",
        ylabel = "Temperature (C)",
        num_points=21
        )
