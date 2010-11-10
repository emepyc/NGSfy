open Getopt
open Types

let model = ref ""
let help = ref false
let compl = ref false
let outfile = ref ""
let infile = ref ""
let prefix = ref ""
let paired = ref false

(* 454 specific *)
let k = ref "0.15"
let lsd = ref "0.1"
let flows = ref "220"

(* Solexa specific *)
let maxerr = ref "0.2"
let maxins = ref "0.0001"
let maxdel = ref "0.0001"

let p_help () = print_endline
"Usage: ngsfy <-m model> <-i infile> <-r repfile> <-o outfile> <-d dictfile> [opts]

Current supported options:
   -m | --model < 454 | Solexa | Illumina >
               The NGS platform to simulate. Solexa and Illumina are synonisms. Required

   -i | --infile < string >
               The input file in fasta format (See below for paired-end options. Required

   -o | --outfile < string >
               The output file.
               The format will be sff or fastq if 454 or Illumina is selected respectively.

   -p | --paired 
               If the reads are paired

   -c | --complement
               If the input sequences are provided in only one strand, this option
               reverse-complement approximately one half of them. Optional (defaults to false)

   -h | --help
               Prints this help and exit

=== 454 options (if -m|--model is set to Solexa|Illumina, these options are ignored ===

   -l | --linker < Titanium | Flx | string >
                The linker in paired-end runs. 

   -k < float >
               Undocummented. Optional (defaults to 0.15)

   -s < float >
               Undocumented. Optional (defaults to 0.1)

   -f | --flows < int >
               Undocummented. Optional (defaults to 220)

=== Solexa/Illumina options (if -m|--model is set to 454, these options are ignored ===

   -e | --maxerror < float >
               Maximum error rates (defaults to 0.2)

   -y | --maxins < float >
               Maximum insertion rates (defaults to 0.0001)

   -z | --maxdel < float >
               Maximum deletion rates (defaults to 0.0001)
"

let specs = 
[
  ('m',"model",None,(atmost_once model (Error "only one mode allowed")));
  ('i',"infile",None,(atmost_once infile (Error "only one input file allowed")));
  ('o',"outfile",None,(atmost_once outfile (Error "only one output file allowed")));
  ('k',"",None, (Some (fun x -> k:=x)));
  ('s',"",None, (Some (fun x -> lsd:=x)));
  ('l',"linker",None,Some (fun x -> match (String.lowercase x) with | "titanium" -> linker := Titanium | "flx" -> linker:=Flx | _ as lk -> linker:=(Custom lk)));
  ('f',"flows",None, Some (fun x -> flows:=x));
  ('c',"complement",(set compl true),None);
  ('h',"help",(set help true), None);
  ('p',"paired",(set paired true),None);
  ('e',"maxerror",None, (Some (fun x -> maxerr:=x)));
  ('y',"maxins",None, (Some (fun x -> maxins:=x)));
  ('z',"maxdels",None, (Some (fun x -> maxdel:=x)));
] 

let main () = 
  parse_cmdline specs print_endline;
  if !help then p_help()
  else
    let io = 
    match String.lowercase !model with
    | "454" -> let o = new Ngsfy.sff_454 !infile !outfile (int_of_string !flows) (float_of_string !k) (float_of_string !lsd) !compl !paired in (o :> Ngsfy.fastq_Solexa)
    | "solexa" | "illumina" -> new Ngsfy.fastq_Solexa !infile !outfile (float_of_string !maxerr) (float_of_string !maxins) (float_of_string !maxdel) !compl !paired
    | _ -> failwith "Error model not supported" 
    in 
    let rec aux () = 
      let v = 
        try Some (io#readRead)
        with
          End_of_file -> None
      in match v with
        Some (h,s) ->
          let ngs_v = io#ngsfyRead s in
          (io#writeRead (h, ngs_v);
           aux())
      | None -> io#close
    in aux ()

let _ = Printexc.print main ()
