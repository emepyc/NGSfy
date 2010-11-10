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

(* Solexa *)
open Types

exception Exit_for of int
let ntIndex = [('A',0);('C',1);('G',2);('T',3);('N',4)]

let errate_Solexa (pos:int):float =
(*  let err = 3e-3 +. 3.3e-8 *. (Gsl_math.pow_int (float_of_int pos) 4 ) in *)
(* alternative for 75bp *)
  let err = 3e-3 +. 1.0e-9 *. (Gsl_math.pow_int (float_of_int pos) 4 ) in
  if err < !maxerror then err else !maxerror

(* let insrate_Solexa = 0.0001 *)
(* let delrate_Solexa = 0.0001 *)

let vec2accum v = 
  try 
    snd (Array.fold_left (
         fun acc x ->
           match acc with
           | ((acu,pos),arr) -> (acu := !acu +. x;
                                 arr.(!pos) <- !acu;
                                 pos := !pos + 1;
                                 acc)
                 ) ((ref 0.,ref 0), Array.make (Array.length v) 0.) v
           )
  with
    Not_found -> failwith "found Not_found"

(* Ordered : A C G T *)
let subs_Solexa = [|
  [|0.000;0.625;0.125;0.250|];
  [|0.470;0.000;0.130;0.400|];
  [|0.100;0.150;0.000;0.750|];
  [|0.200;0.500;0.300;0.000|];
  [|0.0;0.0;0.0;0.0|]
|]

let errQual_Solexa = [|0.001;0.002;0.016;0.0450;0.09;
                        0.130;0.090;0.060;0.050;0.043;
                        0.038;0.036;0.03;0.025;0.0220;
                        0.02;0.019;0.0175;0.015;0.014;
                        0.013;0.012;0.011;0.009;0.009;
                        0.008;0.008;0.007;0.006;0.006;
                        0.005;0.004;0.004;0.004;0.003;
                        0.003;0.003;0.003;0.002;0.002;
                        0.002;0.002;0.002;0.002;0.002;0.1045|]
let cerrQual_Solexa = vec2accum errQual_Solexa

let okQual_Solexa = [|0.000;0.000;0.001;0.012;0.012;
                      0.012;0.012;0.012;0.011;0.011;
                      0.011;0.011;0.011;0.011;0.011;
                      0.011;0.011;0.011;0.011;0.011;
                      0.011;0.011;0.011;0.011;0.011;
                      0.011;0.011;0.011;0.011;0.011;
                      0.010;0.010;0.010;0.010;0.010;
                      0.009;0.009;0.009;0.008;0.008;
                      0.004;0.004;0.004;0.004;0.004;0.584|]
let cokQual_Solexa = vec2accum okQual_Solexa

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

let f_okQual_Solexa() =
  let rval = Gsl_rng.uniform rgen in
  (getRange rval cokQual_Solexa - 5)

let f_errQual_Solexa() =
  let rval = Gsl_rng.uniform rgen in
  (getRange rval cerrQual_Solexa - 5)

let mismatch seq subsvec pos =
  try
  seq.[pos] <- (let randn = Gsl_rng.uniform rgen in
  let accumvec = vec2accum subsvec in  (* Is the compiler smart enough? *)
  match randn with
  | x when x <= accumvec.(0) -> 'A'
  | x when x <= accumvec.(1) -> 'T'
  | x when x <= accumvec.(2) -> 'G'
  | _ -> 'C'
	);
    seq
    with Invalid_argument x -> raise (Invalid_argument "mismatch")
  | Not_found -> raise (Invalid_argument "Not found -- mismatch")

(* Capped at 4 - 34 *)
(* See: http://seqanswers.com/forums/showthread.php?t=4721 *)
let ascii_Solexa1 x = (*print_int x;*)
  if x > 33 then 
(*    (Printf.printf "=> 98 %c\n" (Char.chr 98);Char.chr 98) *)
    Char.chr 98
  else if x < 4 then 
(*    (Printf.printf " => 68 %c\n" (Char.chr 68); Char.chr 68) *)
   Char.chr 68
(*  else (Printf.printf " => %i %c\n" (x+64) (Char.chr (x+64)); Char.chr (x+64)) *)
      else Char.chr (x+64)

let getMismatch_Solexa s q p =
  let newS = mismatch s (subs_Solexa.(List.assoc s.[p] ntIndex)) p in
  try
  begin
    q.[p] <- ascii_Solexa1 (f_errQual_Solexa());
    (newS,q)
  end
  with
    Invalid_argument x -> raise (Invalid_argument "getMismatch_Solexa")
  | Not_found -> raise (Invalid_argument "Not found -- getMismatch_Solexa")


let insert seq pos =
  let newStr = String.create ((String.length seq)+1) in
  begin
    try
    String.blit seq 0 newStr 0 pos;
    String.blit seq pos newStr (pos+1) ((String.length seq) - pos);
    newStr
    with Invalid_argument x -> raise (Invalid_argument "insert")
  end

let getInsert_Solexa s q p =
  let insS = insert s p in
  let newS = mismatch insS (Array.make 4 0.25) p in
  let newQ = insert q p in
  try
  begin
    newQ.[p] <- ascii_Solexa1 (f_errQual_Solexa());
    newQ.[p+1] <- ascii_Solexa1 (f_errQual_Solexa());
(*     print_endline ("\nS -- before:"^s); *)
(*     print_endline ("S -- after :"^newS); *)
(*     print_endline ("Q -- before:"^q); *)
(*     print_endline ("Q -- after :"^newQ); *)
    (newS,newQ)
  end
  with Invalid_argument x -> raise (Invalid_argument "getInsert_Solexa")

let delete seq pos =
  let newStr = String.create ((String.length seq)-1) in
  begin
    try
      String.blit seq 0 newStr 0 pos;
      String.blit seq (pos+1) newStr (pos) ((String.length seq) - pos - 1);
      newStr
    with | Invalid_argument x -> raise (Invalid_argument "delete")
  end

let getDelete_Solexa s q p =
  (delete s p,delete q p)

let applySolexa seq =
   let rec aux acc qual pos =
    assert ((String.length acc) = (String.length qual));
    if pos >= (String.length acc) then FASTQ(acc,qual)
    else
      let perror = errate_Solexa (pos - 1) in 
      let verr = Gsl_rng.uniform rgen in
      let (newS,newQ,newP) = 
        match verr with
        | xerr when xerr <= perror -> let (newS,newQ) = try getMismatch_Solexa acc qual pos with | Invalid_argument x -> raise(Invalid_argument "fst") | Not_found -> raise(Invalid_argument "kkk")  in (newS,newQ,(pos+1))
        | xerr when xerr <= (perror +. !insrate_Solexa) -> let (newS,newQ) = try getInsert_Solexa acc qual pos with Invalid_argument x -> raise(Invalid_argument "snd") in (newS,newQ,(pos+2))
        | xerr when xerr <= (perror +. !insrate_Solexa +. !delrate_Solexa) -> let (newS,newQ) = try getDelete_Solexa acc qual pos with Invalid_argument x -> raise (Invalid_argument "thd") in (newS,newQ,pos)
        | _ -> try (acc,(qual.[pos] <- ascii_Solexa1(f_okQual_Solexa());qual),(pos+1)) with Invalid_argument x -> raise (Invalid_argument "fth")
      in aux newS newQ newP
  in aux seq (String.create (String.length seq)) 0 

(* Compile *)
(* ocamlfind ocamlopt -c simSolexa.ml -package gsl *)
