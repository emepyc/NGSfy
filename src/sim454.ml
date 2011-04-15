(************************************************************************)
(* 
   This module is part of the NGSfy software. NGSfy is free software;
   you can redistribute it and/or modify it under the terms of the
   GNU General Public License as published by the Free Software Foundation;
   either version 2 of the License, or (at your option) any later version.

   This program is distributed as is in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of 
   MERCHANTABILITY of FITNESS FOR A PARTICULAR PURPOSE. See the GNU
   General Public License for more details.

    You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
   USA
*)

(***********************************************************************)
open Types

(* 454 *)
exception No_solution
exception Exit_for of int

(*let k = ref 0.15*) (* User defined -- but defaults to 0.15 -- now in types.ml *)
let k' = 0.15
let linmu1' = 0.2
let linsd1' = 0.1

(* Random number generator *)
let rgen = 
  let rg = Gsl_rng.make (Gsl_rng.default()) in
  begin
    Gsl_rng.set rg (Random.self_init();Random.nativeint (Nativeint.of_int max_int));
    rg
  end

let getRange x vec =
  let prev = ref 0. in
  try (Array.iteri (fun p a -> if x >= !prev & x < a then raise (Exit_for p) else prev:=x) vec;Array.length vec)
  with
  | Exit_for fst -> fst

(* Cutpoint between the first curve (lognormal => negative flow) and the first gaussian (homopolymer length 1 *)
let lognVSgauss () =
  let mu1 = (log linmu1') -. (0.5 *. log(1.+.((linsd1' *. linsd1')/.(linmu1' *. linmu1'))))  in
  let mu2 = 1. in
  let sd1 = sqrt(log(1.+. ((linsd1' *. linsd1') /. (linmu1' *. linmu1')))) in
  let sd2 = k' in
  let rec aux (x1,x2) =
    let x = (x1 +. x2) /. 2. in
    let logn = (1./.(x *. sd1)) *. exp (-.(Gsl_math.pow_int((log x) -. mu1) 2) /. (2. *. sd1 *. sd1)) in
    let gauss = (1./.sd2) *. exp (-.(Gsl_math.pow_int (x -. mu2) 2) /. (2. *. sd2 *. sd2)) in
    let s = gauss -. logn in
    if s < 0.0000001 & s > -.0.0000001 then x
    else if s < 0. then aux (x,x2) else aux (x1,x)
    in aux (linmu1',mu2)

(* Cutpoints between 2 gaussians distribution *)
let gaussVSgauss m1 m2 =
  let mu1 = float_of_int m1 in
  let mu2 = float_of_int m2 in
  let sd1 = k' *. (sqrt mu1) in
  let sd2 = k' *. (sqrt mu2) in
  let a = (sd1 *. sd1 -. sd2 *. sd2) in
  let b = (2. *. mu1 *. (sd2 *. sd2)) -. (2. *. mu2 *. (sd1 *. sd1)) in
  let c = (mu2 *. mu2 *. sd1 *. sd1) -. (mu1 *. mu1 *. sd2 *. sd2) -. (2. *. sd1 *. sd1 *. sd2 *. sd2 *. log(sd1/.sd2)) in
  let sol = Gsl_poly.solve_quadratic a b c in
  match sol with
  | Gsl_poly.Quad_0 -> raise No_solution
  | Gsl_poly.Quad_2 (x1,x2) -> if x1 > mu1 & x2 < mu2 then x1 else x2

let cutpoints = Array.init 30 (
  fun x ->
    if x = 0 then lognVSgauss()
    else
      try gaussVSgauss x (x+1)
      with
      | No_solution -> 0.)

(* let cutpoints = [|0.55958690643310538; 1.4251988012643; 2.46063767946331291; *)
(*     3.47529483059523914; 4.48334859207840442; 5.48844851035663162; *)
(*     6.49196983915027; 7.49454798868530769; 8.49651743467659415; *)
(*     9.49807112229674466; 10.4993282020811272; 11.500366245892323; *)
(*     12.5012379348580556; 13.501980304606036; 14.5026201535369648; *)
(*     15.5031773522121039; 16.5036669466369368; 17.5041005389962478; *)
(*     18.5044872197982; 19.5048342127538561; 20.5051473306435774; *)
(*     21.5054313038014726; 22.5056900209007331; 23.5059267081960215; *)
(*     24.5061440648262767; 25.5063443662526446; 26.5065295442554891; *)
(*     27.5067012494656; 28.5068609007192322; 29.5070097243697305|] *)

let homop s = 
  let lastv = ref (1,s.[0]) in
  let l = String.length s in
  let rec aux p acc =
    if p >= l then List.rev(!lastv::acc)
    else
      if s.[p] = snd(!lastv) && s.[p] != 'N' && s.[p] != 'n' then
        begin
          lastv := ((fst !lastv)+1,(snd !lastv));
          aux (p+1) acc
        end
      else
        let oldv = !lastv in
        begin
          lastv := (1,s.[p]);
          aux (p+1) (oldv::acc)
        end
  in aux 1 []

let expand_i l =
  let rec aux acc l' =
    match l' with
    | (0,_)::tl -> aux acc tl
    | (n,v)::tl -> aux (v::acc) (((n-1),v)::tl)
    | [] -> acc
  in aux [] l

let randomVal r =
  if r = 0 then
    let mu = (log linmu1') -. (0.5 *. log(1.+.((!linsd1 *. !linsd1)/.(linmu1' *. linmu1'))))  in
    let sd = sqrt(log(1. +. ((!linsd1 *. !linsd1) /. (linmu1' *. linmu1')))) in
    Gsl_randist.lognormal rgen mu sd
  else
    let mean = (float_of_int r) in
    let std  = !k *. (sqrt mean) in
    (Gsl_randist.gaussian rgen std) +. mean

let gaussian x mu sd = 1. /. (sd *. (sqrt 2. *. Gsl_math.pi)) *. exp ((-.(Gsl_math.pow_int (x -. mu) 2))/.(2. *. (sd*.sd)))

let lognormal x mu sd = 
  let mu' = (log mu) -. (0.5 *. log(1.+.((sd *. sd)/.(mu *. mu))))  in
  let sd' = sqrt(log(1.+. ((sd *. sd) /. (mu *. mu)))) in
  (1./.(x *. sd')) *. exp (-.(Gsl_math.pow_int((log x) -. mu') 2) /. (2. *. sd' *. sd'))

let probSN x n =
  if n = 0 then
    lognormal x 0.2 0.1
  else
    let mu = (float_of_int n) in
    let sd = k' *. sqrt mu in
    gaussian x mu sd

let probN n = Gsl_math.pow_int 0.25 n

let probNS x l =
  let psn = probSN x l in
  let pn  = probN l in
  let ps = 
    let rec aux n acc =
      if n > (l + 5) then acc   (* Arbitrary limit... max homopolymer length *)
      else
        let psn' = probSN x n in
        let pn'  = probN n in
        aux (n+1) (acc +. (psn' *. pn'))
    in aux 1 (probSN x 0)
  in 
  psn *. pn /. ps

let prob2phred p =
  let phred = int_of_float (floor (-.10. *. (log10 (1. -. p)))) in
    if phred > 40 then 40
    else phred

let appseq2 ctrl slt =
  let vals = Array.init (String.length ctrl) (fun x -> (randomVal 0)/.10.) in
  let rec aux cpos seqlst quals indxs s shifted =
    match seqlst with
    | (n,v)::tl ->
        begin
          try
            let rv = randomVal (if v=ctrl.[cpos] then n else 0) in
            let ppoint = getRange rv cutpoints in
            let p = probNS rv ppoint in
            let q = prob2phred p in
(*             (print_endline ( (string_of_int n) ^ (String.make 1 v) ^ " => " ^ (string_of_float rv) ^ (String.make 1 ctrl.[cpos]) ^ " : " ^ (string_of_int ppoint) ^ " " ^ (string_of_int q)); *)
            vals.(cpos) <- rv;
	    (* ( print_char '\n';print_string "v           : ";print_char v; print_char '\n'; *)
	    (*   print_string "n           : ";print_int n; print_char '\n'; *)
	    (*   print_string "ppoint      : ";print_int ppoint; print_char '\n'; *)
	    (*   print_string "ctrl.[cpos] : ";print_char ctrl.[cpos]; print_char '\n'; *)
	    (*   print_string "rv          : ";print_float rv; print_char '\n'; *)
	    (*   print_string "p           : ";print_float p; print_char '\n'; *)
	    (*   print_string "q           : ";print_int q; print_char '\n'; *)
	    (*   print_string "shifted     : ";print_int shifted; print_char '\n'; *)
            match (ppoint,q) with
	      | (n,0) when n > 0 -> aux (cpos+1) tl ((n,10)::quals) ((n-1,0)::(1,(shifted+1))::indxs) ((n,v)::s) 0 
              | (0,_)
              | (_,0) -> (* What happens if ppoint>1 and q=0?? Here, only one 'N' will be inserted! *)
                if shifted = 2 then
		  aux (cpos+1) tl ((1,0)::quals) ((1,(shifted+1))::indxs) ((1,'N')::s) 0
                else aux (cpos+1) seqlst quals indxs s (shifted+1)
            | (cutp,_) ->
                if cutp != n || ctrl.[cpos] != v then print_endline ("ERR: " ^ (string_of_int n) ^ (String.make 1 v) ^ "x" ^ (string_of_int cutp) ^ (String.make 1 ctrl.[cpos]) ) else ();
(*                aux (cpos+1) tl ((cutp,q)::quals) ((cutp,(shifted+1))::indxs) ((cutp,v)::s) 0 *)
(*              (print_string "INSERTING: "; print_char v; print_char '\n'; *)
		 aux (cpos+1) tl ((cutp,q)::quals) ((cutp-1,0)::(1,(shifted+1))::indxs) ((cutp,v)::s) 0 (*)) (* Lists grow in reverse *)*)
          with
          | Invalid_argument m ->
              print_endline ("Clipping sequence");
              ((List.fold_left (fun acc x -> acc^(String.make 1 x)) "" (expand_i s)),(expand_i quals),(expand_i indxs),vals)
        end
    | [] -> ((List.fold_left (fun acc x -> acc^(String.make 1 x)) "" (expand_i s)),(expand_i quals),(expand_i indxs),vals)
  in aux 0 slt [] [] [] 0

let apply454 s ctrlseq =
    let slst = homop s in
    let (newseq,quals,indxs,vals) = appseq2 ctrlseq slst in
    begin
      let vals = Array.to_list vals in
      assert (List.length indxs = String.length newseq);
      assert (List.length quals = String.length newseq);
      assert (List.length vals = String.length ctrlseq);
      SFF454 { sequence=newseq;
                signals=vals;
                indxs=indxs;
                quals=quals
              }
    end


(* To compile *)
(* ocamlfind ocamlopt -c sim454.ml -package gsl *)
