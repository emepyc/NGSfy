external cleanseq : string -> string = "cleanseq"
external ltrim : string -> string = "ltrim"
external revcomp : string -> string = "revcomp"
external cleanseq_light : string -> string = "cleanseq_light"
val levdistance : string -> string -> int
val wordSimilarity : string -> string -> float
(* misspelled returns if both strings seems to be the same string (allowing misspelleds) *)
val sameString : string -> string -> bool
(* str2wordList : CAUTION! This function uses Genlex and Stream, and fails with things like 2'*)
(* or in general, if a single quote is used at the beginning of a token *)
val str2wordList : string -> string list
(* splitInBlanks *)
val splitInBlanks : string -> string list
