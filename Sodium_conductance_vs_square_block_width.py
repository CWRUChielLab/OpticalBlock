import nrnaxon
from neuron import h

def sodium_block_conductance_at_width(block_width):
    axon_length = 300.
    axon_num_sections = 3000

    def sodium_block(gNa):
        axon = nrnaxon.Axon(length=axon_length, num_sections=axon_num_sections)
        axon.insert_stim()
        axon.set_gNa(gNa,
                (axon_length - block_width)/2.,
                (axon_length + block_width)/2.
                )
        return nrnaxon.is_blocked(axon)

    bounds = nrnaxon.bisect(sodium_block, 0., 8e-3, 20)
    return sum(bounds)/2.


if __name__ == "__main__":
    nrnaxon.record_plot(
        func=sodium_block_conductance_at_width,
        xmin = 42,
        xmax = 82,
        title = "Sodium conductance vs square block width",
        xlabel = "Block width (um)",
        ylabel = "Sodium conductance (S per square cm)",
        num_points=21
        )
