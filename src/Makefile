VERSION = 0.01

# Ocaml distribution packages needed
REQUIRES = extlib str getopt gsl

# Own libraries included
INCLPATH = .
INCLUDES = $(INCLPATH)/myString

# Common options used
COMMONOPTS = -ccopt -static -package "$(REQUIRES)" -I $(INCLUDES)

# Compile with ocamlfind
OCAMLFIND = ocamlfind
OCAMLOPT = $(OCAMLFIND) ocamlopt $(COMMONOPTS)
OCAMLDEP = $(OCAMLFIND) ocamldep $(COMMONOPTS)
OCAMLLEX = ocamllex
CCLIBS = $(INCLUDES)/studlibs


OBJS = \
	bin_io.cmx \
	types.cmx \
	fasta.cmx \
	sff.cmx \
	sim454.cmx \
	simSolexa.cmx \
	ngsfy.cmx \



all: ngsfy

bin_io.cmx: bin_io.ml
	$(OCAMLOPT) -c $<

types.cmx: types.ml 
	$(OCAMLOPT) -c $<

fasta.ml: fasta.mll
	$(OCAMLLEX) $<

fasta.cmx: fasta.ml
	$(OCAMLOPT) -I $(INCLUDES) -c $<

sff.cmx: sff.ml
	$(OCAMLOPT) -c $<

sim454.cmx: sim454.ml
	$(OCAMLOPT) -c $<

simSolexa.cmx: simSolexa.ml
	$(OCAMLOPT) -c $<

ngsfy.cmx: ngsfy.ml
	$(OCAMLOPT) -I $(INCLUDES) -c $<

NGSfy.cmx: NGSfy.ml
	$(OCAMLOPT) -c $<

ngsfy: $(OBJS) NGSfy.cmx
	$(OCAMLOPT) -cclib $(CCLIBS)/lib_wrap_stubs.a -linkpkg $(OBJS) NGSfy.cmx -o ngsfy


clean:
	rm -f *.cm[aiox] *.o *.so *.a *.cmxa
	rm -f fasta.ml
