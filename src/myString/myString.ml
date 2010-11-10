external cleanseq : string -> string = "cleanseq"
external cleanseq_light : string -> string = "cleanseq_light"
external ltrim : string -> string = "ltrim"
external revcomp : string -> string = "revcomp"

open ExtLib

(* http://en.wikibooks.org/wiki/Algorithm_implementation/Strings/Levenshtein_distance#OCaml *)

let explode w =
	let rec aux pos lst =
		if pos = String.length w then lst
    else aux (succ pos) (w.[pos]::lst)
  in List.rev (aux 0 [])
;;

(* Minimum of three integers. This function is deliberately
 * not polymorphic because (1) we only need to compare integers 
 * and (2) the OCaml compilers do not perform type specialization 
 * for user-defined functions. *)
let rec minimum (x : int) y z =
  if y < x then minimum y x z else
  if z < y then minimum x z y else x
 
(* Matrix initialization. *)
let init_matrix n m =
  let init_col = Array.init m in
  Array.init n (function
    | 0 -> init_col (function j -> j)
    | i -> init_col (function 0 -> i | _ -> 0)
  )
 
(* Computes the Levenshtein distance between two unistring.
 * If you want to run it faster, add the -unsafe option when
 * compiling or use Array.unsafe_* functions (but be carefull 
 * with these well-named unsafe features). *)
let distance_utf8 x y =
  let n = Array.length x and m = Array.length y in
  match n, m with
  | 0, _ -> m
  | _, 0 -> n
  | _ -> let matrix = init_matrix (n + 1) (m + 1) in
    for i = 1 to n do
      let s = matrix.(i) and t = matrix.(i - 1) in
      for j = 1 to m do
        let cost = abs (compare x.(i - 1) y.(j - 1)) in
        s.(j) <- minimum (t.(j) + 1) (s.(j - 1) + 1) (t.(j - 1) + cost)
      done
    done;
    matrix.(n).(m)
 
(* This function takes two strings, convert them to unistring (int array)
 * and then call distance_utf8, so we can compare utf8 strings. Please
 * note that you need Glib (see LablGTK). *)
let levdistance x y =
  distance_utf8 (Array.of_list (explode x)) (Array.of_list (explode y))

let wlexer = Genlex.make_lexer [",";".";";";"[";"]";"\"";"(";")";"'"]

(* CAUTION! This function uses Genlex and Stream, and fails with things like 2' *)
(* or in general, if a single quote is used at the beginning of a token *)
let str2wordList str =
	(print_endline str;
	let stream = wlexer (Stream.of_string str) in
	let rec aux acc =
		let nextWord =
			try
				Some (Stream.next stream)
			with
				| Stream.Failure -> None
		in match nextWord with
			| Some (Genlex.Ident x) -> aux (x::acc)
			| Some (Genlex.Int i) -> aux ((string_of_int i)::acc)
			| Some (Genlex.Float f) -> aux ((string_of_float f)::acc)
			| None -> List.sort ~cmp:Pervasives.compare acc
			| _ -> aux acc
	in aux []
)

let splitInBlanks str = String.nsplit str " "

(* Cambiar para que admita UTF-8 *)
let wordSimilarity str1 str2 =
	let (wl1,wl2) = (str2wordList str1, str2wordList str2) in
	let rec aux l1 l2 acc =
		match (l1,l2) with
			| (w1::t1,w2::t2) when w1 = w2 -> aux t1 t2 (succ acc)
			| (w1::t1,w2::t2) -> max (aux l1 t2 acc) (aux t1 l2 acc)
			| (_,_) -> (float_of_int acc) /. float_of_int ((max (List.length wl1) (List.length wl2)))
	in aux wl1 wl2 0

(* The following functios try to ditinguish if 2 strings are the same -- allowing common typos *)
let hasNumbers w =
	let rec aux l =
		match l with
			| [] -> false
			| x::xs -> let xcode = Char.code x in if xcode >= 48 && xcode <= 57 then true
								 else aux xs
	in aux w

let isSubString wl1 wl2 =
	let swl1 = List.sort ~cmp:Pervasives.compare wl1
	and swl2 = List.sort ~cmp:Pervasives.compare wl2 in
	let rec aux x y acc =
		match x,y with
			| [],[] -> acc
			| w,[]
			| [],w -> w@acc
			| c1::c1s, c2::c2s -> 
				match (compare c1 c2) with
					| -1 -> aux c1s y (c1::acc)
					| 0 -> aux c1s c2s acc
					| 1 | _ -> aux x c2s (c2::acc)
	in aux swl1 swl2 []
	

let sameString : string -> string -> bool =
	fun s1 s2 -> 
		if s1 = s2 then true else
		if abs ((String.length s1) - (String.length s2)) > 4 then
			false
		else
(*			let (wl1,wl2) = (str2wordList s1, str2wordList s2) in *) (* See note on function str2wordList definition *)
			let (wl1,wl2) = (splitInBlanks s1, splitInBlanks s2) in
			let diffWLength = abs ((List.length wl1) - (List.length wl2)) in
			if diffWLength > 1 then
				false
			else
				let diffLetters = isSubString (List.concat (List.map (fun x -> explode x) wl1)) (List.concat (List.map (fun x -> explode x) wl2)) in
				if List.length diffLetters > 3 then false
				else
					if (hasNumbers diffLetters) then false else true


