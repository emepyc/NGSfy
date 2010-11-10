

open Types
open Bin_io

let titanium_linker = "TCGTATAACTTCGTATAATGTATGCTATACGAAGTTATTACG"
let flx_linker = "GTTGGAACCGAAAGGGTTTGAATTCAAACCCTTTCGGTTCCAAC" (* palindromic *)

(* let parse_header header = *)
(*   let tokenfirst = Str.split (Str.regexp " ") header in *)
(*   match tokenfirst with *)
(*   | [] -> failwith "Can't parse header" *)
(*   | fs::rest -> let tokensL = Str.split (Str.regexp "|") fs in *)
(*     let tokensA = Array.of_list tokensL in *)
(*     { abbrev = tokensA.(0); *)
(*       chrom = tokensA.(1); *)
(*       from = int_of_string tokensA.(2); *)
(*       length = int_of_string tokensA.(3); *)
(*       str = '+' *)
(*     }  *)

(*let rinfo2key (r:r_info):string = r.abbrev ^ r.chrom ^ (string_of_int r.from) ^ (string_of_int r.length) *)

(* let procL (s:string):(r_info * r_info) = *)
(*   let ffs = Str.split (Str.regexp "\t") s in *)
(*   let affs = Array.of_list ffs in *)
(* (\*  let h = parse_header affs.(0) in *\) *)
(*   let h = affs.(0) in *)
(*   let fis = Str.split (Str.regexp "|") affs.(1) in *)
(*   let afis = Array.of_list fis in *)
(*   let (b',e') = (int_of_string affs.(2),int_of_string affs.(3)) in *)
(* (\*  let (b,l,str) = if b' < e' then (b',e'-b'+1,'+') else (e',b'-e'+1,'-') in *\) *)
(*   (\* We correct an error in reps files that gives n-1 bases in the coords *\) *)
(*   let (b,l,str) = if b' < e' then (b',e'-b' ,'+') else (e',b'-e' ,'-') in *)
(*   let rinfo = { abbrev = afis.(1); *)
(*                 chrom = afis.(0); *)
(*                 from = b; *)
(*                 length = l; *)
(*                 str = str *)
(*               } in *)
(*   (h,rinfo) *)

(* let load_reps (ifn:string):(string,r_info) Hashtbl.t =  *)
(*   let htbl = Hashtbl.create 1000000 in *)
(*   let ifh = open_in ifn in *)
(*   let rec aux () =  *)
(*     let line =  *)
(*       try *)
(*         Some (input_line ifh) *)
(*       with *)
(*       | End_of_file -> None *)
(*     in *)
(*     match line with *)
(*     | None -> close_in ifh; htbl *)
(*     | Some l -> let (k,v) = procL l in aux (Hashtbl.add htbl (rinfo2key k) v); *)
(*   in aux () *)


(* let header_generator(():unit):(unit->string) = *)
(*   let v = ref 0 in *)
(*   fun () -> let n = !v in *)
(*   ( v := !v + 1; *)
(*    Printf.sprintf "%09d" n ) *)

class virtual platform f_in comp paired =
  let lexer = Lexing.from_channel (open_in f_in) in
(*  let d = new Headers.odict dict_name in *)
  object (self)
(*    val mutable f = header_generator() *)
(*    val dict = d *)
    val compl = comp
(*    val reps = load_reps repfile *)
    val paired = paired
(*    method private nextName = f() *)
    method virtual ngsfyRead : seqType -> read
    method virtual writeRead : (headType * read) -> unit
    method virtual close : unit
    method make_se h s do_rc =
(*      let rinfo = parse_header h in  *)
      let s_clean = MyString.cleanseq s in
      (SEH h,SES (if do_rc then MyString.revcomp s_clean else s_clean))
    method virtual make_pe : (string*string) -> (string*string) -> bool -> (headType*seqType)
    method readRead = (* Lexing.lexbuf -> headType*seqType *)
      let (h,s) = Fasta.pfasta lexer in
      let do_comp = if compl & ((Random.float 1.) < 0.5) then true else false in
      try 
        if paired then (* File is paired - end*)
          let (h',s') = Fasta.pfasta lexer in
          self#make_pe (h,h') (s,s') do_comp
        else
          self#make_se h s do_comp
      with Failure x -> raise (Failure "caught")
    initializer
      Random.self_init();
  end

class fastq_Solexa f_in f_out maxerr maxins maxdel comp paired =
  object (self)
    inherit platform f_in comp paired
    val cout = open_out f_out
    method make_pe (h1,h2) (s1,s2) do_rc = 
(*      let r1info = parse_header h1 in *)
(*      let r2info = parse_header h2 in *)
      (PEH(h1,h2),PESol((if do_rc then MyString.revcomp s1 else s1) ,(if do_rc then s2 else MyString.revcomp s2)))
    method ngsfyRead k = (* seqtype -> read *)
      match k with
      | SES s -> SimSolexa.applySolexa s
      | PESol (s1,s2) ->
          let ns1 = SimSolexa.applySolexa s1 in
          let ns2 = SimSolexa.applySolexa s2 in
          FASTQpe (ns1,ns2)
      | _ -> failwith "Error - Wrong seq provided"
    method writeRead (h,sq) =
      let hnude = match h with | SEH (hh) -> hh | _ -> failwith "Error right here!" in
(*      let codeH = prefix ^ self#nextName in *)
      match h,sq with
      | (SEH x),(FASTQ(s,q)) -> 
(*          let rs = Hashtbl.find_all reps (rinfo2key rinfo) in *)
(*          dict#writeDict (rinfo::rs); *)
          Printf.fprintf cout "@%s\n%s\n+\n%s\n" x s q
      | (PEH(_,_)),(FASTQpe(FASTQ(s1,q1),FASTQ(s2,q2))) ->
(*          ignore(self#nextName); *)
(*          let rs1 = Hashtbl.find_all reps (rinfo2key r1info) in *)
(*          let rs2 = Hashtbl.find_all reps (rinfo2key r2info) in *)
(*          dict#writeDict (r1info::rs1); *)
(*          dict#writeDict (r2info::rs2); *)
          Printf.fprintf cout "@%s\n%s\n+\n%s\n" (hnude^"/1") s1 q1;
          Printf.fprintf cout "@%s\n%s\n+\n%s\n" (hnude^"/2") s2 q2
      | _ -> failwith "writeRead(Solexa) -- Invalid type of header or sequence"
    method close = Pervasives.close_out cout
    initializer
      maxerror := maxerr;
      insrate_Solexa := maxins;
      delrate_Solexa := maxdel;
  end
    
class sff_454 f_in f_out rounds vk vlsd comp paired =
  object (self)
    inherit platform f_in comp paired as super'
    inherit Sff.outsff f_out rounds as super
    val rounds = rounds
    method readRead = 
      match super'#readRead with
      | (_ as e,PE454 s) -> (e,PE454 (key^s))
      | (_ as e,SES s) -> (e,SES (key^s))
      | _ -> failwith "Impossible in readRead(454)"
    method make_pe (h1,h2) (s1,s2) do_rc =
      let lkr = 
        match !linker with
        | NoLink -> failwith "No linker provided for paired-ends reads, but paired-end data found in input file"
        | Titanium -> titanium_linker
        | Flx -> flx_linker
        | Custom x -> x
      in
      let linker_l = String.length lkr in
      let seq_l = String.length s1 in
(*      let s2_rc = MyString.revcomp s2 in *) (* |--- s2 -->|+++linker+++|--- s1 --->| *)
      let newS = s2 ^ lkr ^ s1 in (* |--- s2 -->|+++linker+++|--- s1 --->| *)
      let b = (Random.int (seq_l - linker_l)) + linker_l in
      try
        let newSS = String.sub newS b seq_l in
(*        let r1info = parse_header h1 in *)
(*        let r2info = parse_header h2 in *)
        (PEH(h1,h2),PE454 (if do_rc then MyString.revcomp newSS else newSS))
      with
      | Invalid_argument "String.sub" -> raise (Invalid_argument "caught2")
      
    method ngsfyRead k = 
      match k with
      | SES s -> Sim454.apply454 s ctrlseq
      | PE454 s -> Sim454.apply454 s ctrlseq
      | _ -> failwith "Error - Wrong seq provided"
    method writeRead (h,read) = 
(*      let codeH = prefix ^ self#nextName in *)
      let hnude = match h with | SEH(hh) -> hh | _ -> failwith "Error here again!" in
      match (h,read) with
      | (SEH h), (SFF454 readR) ->
(*          let rs = Hashtbl.find_all reps (rinfo2key rinfo) in  *)
(*          dict#writeDict (rinfo::rs); *)
          super#writeRead_h (h,read)
      | (PEH(h1,h2),SFF454 readR) ->
(*          ignore(self#nextName); *)
(*          let rs1 = Hashtbl.find_all reps (rinfo2key r1info) in *)
(*          let rs2 = Hashtbl.find_all reps (rinfo2key r2info) in *)
(*          dict#writeDict (r1info::rs1); *)
(*          dict#writeDict (r2info::rs2); *)
          super#writeRead_h (hnude,read)
      | _ -> failwith "writeRead(454) -- Invalid type of header or sequence"
    initializer
      k := vk;
      linsd1 := vlsd;
  end

(********)
(* let sffio = new sff_454 "/home/pignatelli/tmp/FW1F50V01.sff" "/home/pignatelli/tmp/yyy.sff";; *)
(* let v = sffio#readRead;; *)
(* sffio#writeRead v;; *)
(* sffio#close_out;; *)
