OUTPUTS = demo_block.pdf

all : $(OUTPUTS)

clean:
	rm -f $(OUTPUTS) $(OUTPUTS:%.pdf=%.csv)

%.csv : %.yaml nrnaxon.py
	nrniv -python nrnaxon.py $<

%.pdf : %.csv
	./plot.R $< $@
