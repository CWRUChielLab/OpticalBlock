OUTPUTS = demo_block.pdf
CSVFILES = $(OUTPUTS:%.pdf=%.csv)
.PRECIOUS: $(CSVFILES)

all : $(OUTPUTS)

clean:
	rm -f $(OUTPUTS) $(CSVFILES)

%.csv : %.yaml nrnaxon.py
	nrniv -python nrnaxon.py $<

%.pdf : %.yaml %.csv
	./plot.R $^ $@
