VERSION = 0.01
REQUIRES = extlib
COMMONOPTS = -package "$(REQUIRES)"
OCAMLFIND = ocamlfind
OCAMLOPT = $(OCAMLFIND) ocamlopt $(COMMONOPTS)
OCAMLDEP = $(OCAMLFIND) ocamldep $(COMMONOPTS)
OCAMLMKLIB = ocamlmklib -custom
OCAMLC= ocamlc.opt

native: myString.cmxa
bytecode: myString.cma

strings_stub.o: strings_stub.c
	$(OCAMLC) -c -verbose $<

lib_wrap_stubs.a: strings_stub.o
	$(OCAMLMKLIB) -o _wrap_stubs  $<
	mv lib_wrap_stubs.a studlibs/

myString.cmi: myString.mli
	$(OCAMLOPT) -c $<

myString.cmx: myString.ml myString.cmi
	$(OCAMLOPT) -c $<

myString.cmxa:  myString.cmx  lib_wrap_stubs.a
	$(OCAMLOPT) -a  -o $@  $<  -cclib -l_wrap_stubs

myString.cmo: myString.ml myString.cmi
	$(OCAMLC) -c $<

myString.cma:  myString.cmo  lib_wrap_stubs.a
	$(OCAMLC) -a  -o $@  $<  -dllib -l_wrap_stubs


distclean: clean
clean:
	rm -f *.cm[aiox] *.o *.so *.a *.cmxa

